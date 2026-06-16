"""Local SRA discovery and pilot-database helpers.

This module supports the prototype question "can we work locally first?" by
finding SRA-derived files on disk and creating small sampled BLAST databases
from existing FASTA or local SRA files.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re
import subprocess

from config import resource_root, runtime_data_dir, tool_name
from database_registry import (
    create_database_from_fasta,
    get_database_by_prefix,
    register_existing_database,
    slugify,
)
from database_size import database_storage_bytes, format_bytes


FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}
SRA_SUFFIX = ".sra"


@dataclass(frozen=True)
class SraFileSummary:
    path: str
    size_bytes: int
    size_label: str


@dataclass(frozen=True)
class SraProject:
    accession: str
    root_path: str
    sra_files: list[SraFileSummary]
    fasta_files: list[SraFileSummary]
    blast_prefixes: list[str]
    blast_database_bytes: int
    blast_database_size_label: str
    total_bytes: int
    total_size_label: str
    status: str


def file_summary(path: Path) -> SraFileSummary:
    """Return path and size metadata for one local file."""
    size_bytes = path.stat().st_size
    return SraFileSummary(
        path=str(path),
        size_bytes=size_bytes,
        size_label=format_bytes(size_bytes),
    )


def configured_sra_roots() -> list[Path]:
    """Return the local folders that the SRA workbench should scan."""
    roots: list[Path] = []
    env_roots = os.environ.get("COBLAST_SRA_DIR") or os.environ.get("SRA_DATA_DIR")
    if env_roots:
        roots.extend(Path(part).expanduser() for part in env_roots.split(os.pathsep) if part)

    roots.append(runtime_data_dir() / "sra")
    roots.append(resource_root().parent / "SRA_data")

    unique: list[Path] = []
    seen: set[str] = set()
    for root in roots:
        try:
            normalized = str(root.expanduser().resolve())
        except OSError:
            normalized = str(root.expanduser())
        if normalized not in seen:
            seen.add(normalized)
            unique.append(Path(normalized))
    return unique


def sra_toolkit_bin() -> Path | None:
    """Find a local SRA Toolkit bin directory, if one is configured or nearby."""
    candidates: list[Path] = []
    env_bin = os.environ.get("SRA_TOOLKIT_BIN")
    if env_bin:
        candidates.append(Path(env_bin).expanduser())

    sibling_root = resource_root().parent / "sratoolkit"
    if sibling_root.exists():
        candidates.extend(path / "bin" for path in sibling_root.glob("sratoolkit.*"))

    for candidate in candidates:
        try:
            resolved = candidate.resolve()
        except OSError:
            resolved = candidate
        if (resolved / tool_name("fastq-dump")).exists():
            return resolved
    return None


def sra_tool_exe(name: str) -> Path:
    """Resolve an SRA Toolkit executable for local pilot conversion."""
    bin_dir = sra_toolkit_bin()
    if bin_dir is None:
        raise FileNotFoundError(
            "Could not find SRA Toolkit. Set SRA_TOOLKIT_BIN to the toolkit bin directory."
        )
    exe_path = bin_dir / tool_name(name)
    if not exe_path.exists():
        raise FileNotFoundError(f"Could not find {exe_path}.")
    return exe_path


def project_accession(path: Path) -> str:
    """Choose a readable accession/project name from an SRA project folder."""
    sra_files = list(path.glob(f"*{SRA_SUFFIX}"))
    if sra_files:
        return sra_files[0].stem
    return path.name


def find_fasta_files(files: list[Path]) -> list[Path]:
    """Select FASTA files from a pre-collected file list."""
    return sorted(path for path in files if path.suffix.lower() in FASTA_SUFFIXES)


def find_sra_files(files: list[Path]) -> list[Path]:
    """Select SRA files from a pre-collected file list."""
    return sorted(path for path in files if path.suffix.lower() == SRA_SUFFIX)


def find_blast_prefixes(files: list[Path]) -> list[str]:
    """Find likely BLAST database prefixes from a pre-collected file list."""
    prefixes: list[Path] = []
    for path in files:
        suffix = path.suffix.lower()
        if suffix in {".nal", ".pal"}:
            prefixes.append(path.with_suffix(""))
        elif suffix in {".nin", ".pin"} and not re.search(r"\.\d+$", path.stem):
            prefixes.append(path.with_suffix(""))

    unique: list[str] = []
    seen: set[str] = set()
    for prefix in prefixes:
        normalized = str(prefix)
        if normalized not in seen:
            seen.add(normalized)
            unique.append(normalized)
    return unique


def summarize_project(root: Path) -> SraProject:
    """Summarize SRA, FASTA, and BLASTDB artifacts under one folder."""
    # One filesystem walk feeds every breakdown below.
    all_files = [path for path in root.rglob("*") if path.is_file()]

    sra_files = find_sra_files(all_files)
    fasta_files = find_fasta_files(all_files)
    blast_prefixes = find_blast_prefixes(all_files)
    blast_database_bytes = sum(database_storage_bytes(prefix) for prefix in blast_prefixes)

    total_bytes = 0
    for path in all_files:
        try:
            total_bytes += path.stat().st_size
        except OSError:
            continue

    if blast_prefixes:
        status = "blast-ready"
    elif fasta_files:
        status = "fasta-ready"
    elif sra_files:
        status = "sra-only"
    else:
        status = "empty"

    return SraProject(
        accession=project_accession(root),
        root_path=str(root),
        sra_files=[file_summary(path) for path in sra_files],
        fasta_files=[file_summary(path) for path in fasta_files],
        blast_prefixes=blast_prefixes,
        blast_database_bytes=blast_database_bytes,
        blast_database_size_label=format_bytes(blast_database_bytes),
        total_bytes=total_bytes,
        total_size_label=format_bytes(total_bytes),
        status=status,
    )


def discover_sra_projects() -> list[SraProject]:
    """Scan configured SRA roots and return project summaries."""
    projects: list[SraProject] = []
    seen_roots: set[str] = set()
    for root in configured_sra_roots():
        if not root.exists():
            continue

        candidate_dirs = [path for path in root.iterdir() if path.is_dir()]
        if any(path.suffix.lower() == SRA_SUFFIX for path in root.iterdir() if path.is_file()):
            candidate_dirs.append(root)

        for candidate in candidate_dirs:
            normalized = str(candidate.resolve())
            if normalized in seen_roots:
                continue
            seen_roots.add(normalized)
            project = summarize_project(candidate)
            if project.sra_files or project.fasta_files or project.blast_prefixes:
                projects.append(project)

    return sorted(projects, key=lambda project: project.accession.lower())


def source_fasta_for_blast_prefix(
    db_prefix_path: str | Path, fasta_files: list[SraFileSummary]
) -> str:
    """Choose the source FASTA that most likely produced a discovered BLAST DB."""
    if not fasta_files:
        return ""

    prefix_name = Path(db_prefix_path).name.casefold()
    same_stem = [
        fasta.path
        for fasta in fasta_files
        if Path(fasta.path).stem.casefold() == prefix_name
    ]
    if len(same_stem) == 1:
        return same_stem[0]
    if len(fasta_files) == 1:
        return fasta_files[0].path
    return ""


def copy_fasta_subset(source_fasta_path: str | Path, output_fasta_path: str | Path, max_records: int) -> int:
    """Copy the first N FASTA records into a smaller pilot FASTA."""
    if max_records < 1:
        raise ValueError("Pilot record count must be at least 1.")

    source = Path(source_fasta_path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"Source FASTA does not exist: {source}")

    output = Path(output_fasta_path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    records_written = 0
    current_record_active = False
    with source.open("r", encoding="utf-8", errors="replace") as reader:
        with output.open("w", encoding="utf-8") as writer:
            for line in reader:
                if line.startswith(">"):
                    if records_written >= max_records:
                        break
                    records_written += 1
                    current_record_active = True
                if current_record_active:
                    writer.write(line)

    if records_written == 0:
        raise ValueError(f"No FASTA records were found in {source}.")
    return records_written


def create_pilot_database_from_fasta(
    *,
    accession: str,
    source_fasta_path: str | Path,
    max_records: int,
):
    """Create and register a sampled nucleotide BLAST database from FASTA."""
    accession_slug = slugify(accession, default="sra_project")
    pilot_dir = runtime_data_dir() / "sra_pilots" / f"{accession_slug}_{max_records}"
    pilot_fasta = pilot_dir / f"{accession_slug}_{max_records}.fasta"
    db_prefix = pilot_dir / "blastdb" / f"{accession_slug}_{max_records}"
    records_written = copy_fasta_subset(source_fasta_path, pilot_fasta, max_records)

    return create_database_from_fasta(
        display_name=f"SRA pilot {accession} ({records_written} reads)",
        db_type="nucl",
        source_fasta_path=pilot_fasta,
        db_prefix_path=db_prefix,
        description=f"Pilot nucleotide BLAST database sampled from {accession}.",
        category="sra",
        notes=(
            f"Pilot database created from the first {records_written} FASTA records. "
            "Use for local runtime simulation before full SRA analysis."
        ),
    )


def convert_sra_to_pilot_fasta(
    *,
    accession: str,
    sra_path: str | Path,
    max_spots: int,
) -> Path:
    """Use fastq-dump to convert a limited number of local SRA spots to FASTA."""
    if max_spots < 1:
        raise ValueError("Pilot spot count must be at least 1.")

    source = Path(sra_path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"SRA file does not exist: {source}")

    accession_slug = slugify(accession or source.stem, default="sra_project")
    output_dir = runtime_data_dir() / "sra_pilots" / f"{accession_slug}_{max_spots}" / "fasta"
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        str(sra_tool_exe("fastq-dump")),
        "--fasta",
        "0",
        "--skip-technical",
        "--readids",
        "--maxSpotId",
        str(max_spots),
        "--outdir",
        str(output_dir),
        str(source),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError((completed.stderr or completed.stdout).strip())

    fasta_files = sorted(path for path in output_dir.iterdir() if path.suffix.lower() in FASTA_SUFFIXES)
    if not fasta_files:
        raise RuntimeError(f"fastq-dump completed but no FASTA file was created in {output_dir}.")
    return fasta_files[0]


def register_sra_blast_database(
    *,
    accession: str,
    db_prefix_path: str | Path,
    source_fasta_path: str | Path | None = None,
):
    """Register an existing SRA-derived nucleotide BLAST database."""
    source = source_fasta_path
    if not source:
        existing = get_database_by_prefix(db_prefix_path)
        if existing is not None and existing.source_fasta_path:
            source = existing.source_fasta_path
    return register_existing_database(
        display_name=f"SRA {accession} reads",
        db_type="nucl",
        db_prefix_path=db_prefix_path,
        source_fasta_path=source,
        description=f"Local nucleotide BLAST database prepared from {accession}.",
        category="sra",
        notes="Registered from the SRA workbench.",
    )
