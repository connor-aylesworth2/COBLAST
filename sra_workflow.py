"""Local SRA discovery and pilot-database helpers.

This module supports the prototype question "can we work locally first?" by
finding SRA-derived files on disk and creating small sampled BLAST databases
from existing FASTA or local SRA files.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

from config import blast_exe, ensure_tool_bin, resource_root, runtime_data_dir, tool_name
from database_registry import (
    create_database_from_fasta,
    get_database_by_prefix,
    register_existing_database,
    slugify,
    verify_database_prefix,
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
    """Resolve an SRA Toolkit executable, auto-installing the toolkit if absent."""
    bin_dir = ensure_tool_bin("sra", sra_toolkit_bin)
    if bin_dir is None:
        raise FileNotFoundError(
            "Could not find SRA Toolkit. Set SRA_TOOLKIT_BIN to the toolkit bin directory."
        )
    exe_path = bin_dir / tool_name(name)
    if not exe_path.exists():
        raise FileNotFoundError(f"Could not find {exe_path}.")
    return exe_path


# Run-level accessions only (SRR/ERR/DRR). Study/experiment/sample accessions
# (SRP/SRX/SRS/PRJ...) are rejected so a tester can't accidentally queue a whole
# multi-terabyte study; prefetch pulls the actual per-sample reads from a run.
RUN_ACCESSION_RE = re.compile(r"^[SED]RR\d+$")


def sra_download_dir() -> Path:
    """The scanned SRA root that new prefetch downloads should land in."""
    target = runtime_data_dir() / "sra"
    target.mkdir(parents=True, exist_ok=True)
    return target


def parse_run_accessions(raw: str) -> list[str]:
    """Split user input into de-duplicated, validated SRA run accessions."""
    tokens = [token.strip().upper() for token in re.split(r"[\s,;]+", raw) if token.strip()]
    if not tokens:
        raise ValueError("Enter at least one SRA run accession, e.g. SRR123456.")
    bad = [token for token in tokens if not RUN_ACCESSION_RE.match(token)]
    if bad:
        raise ValueError(
            "Only SRA run accessions (SRR/ERR/DRR) are allowed. Rejected: "
            + ", ".join(bad)
        )
    seen: set[str] = set()
    unique: list[str] = []
    for token in tokens:
        if token not in seen:
            seen.add(token)
            unique.append(token)
    return unique


def build_fetch_script_lines(
    accessions: list[str],
    sra_dir: Path,
    prefetch: Path,
    fasterq_dump: Path,
    makeblastdb: Path,
) -> list[tuple[str, str]]:
    """Per accession: prefetch the .sra, expand it to FASTA, index it as a blastdb.

    Returns (progress_header, command) pairs. The terminal echoes the header
    before each command so the user always sees "run i/N, step x/3". Steps 1 and 2
    each also show a live native progress bar; step 3 (makeblastdb) has no progress
    option at all, so the header is the only position signal there. The header is
    plain ASCII with no apostrophes so it single-quotes cleanly in PowerShell and sh.

    Three steps so a fetched run lands *blast-ready* rather than stalling at
    fasta-ready: the workbench only shows a "Register BLASTDB" control once
    .nin/.nal index files exist, so without the makeblastdb step a tester ends up
    with a FASTA they cannot register (0 B BLASTDB) and cannot feed to eToL.
      1. prefetch --progress --max-size u: download the .sra with a live progress
         bar and no size cap. prefetch's default cap is 20G and would silently
         truncate a larger run; the accession guard already blocks whole studies,
         so the cap only hurts here.
      2. fasterq-dump --fasta --split-spot --skip-technical --progress: expand the
         whole run to one FASTA with a live progress bar (fastq-dump has none) and
         multi-threaded, so a big run's convert isn't a long silent stall. Every
         biological read is its own record so paired mates stay separate instead of
         concatenating into chimeric sequences. --seq-defline '$ac.$si.$ri'
         reproduces fastq-dump --readids' unique <acc>.<spot>.<read> ids: without a
         read id in the defline, --split-spot gives both mates of a spot the SAME id
         and step 3's -parse_seqids collides/drops them. Single-quoted so neither
         PowerShell nor sh expands the $ tokens. -O writes <acc>.fasta beside the
         .sra; -t keeps fasterq-dump's multi-GB scratch on the same disk (its
         default is the cwd, which the terminal launches from); -f overwrites a
         partial FASTA left by a re-run.
      3. makeblastdb -parse_seqids: index that FASTA into a nucleotide blastdb
         beside it. -parse_seqids builds the id index eToL read recovery and the
         human filter (blastdbcmd) depend on; without it recovery silently degrades.
    Double-quoted paths work in both cmd.exe and POSIX shells; accessions are
    pre-validated (SRR/ERR/DRR) so they never need quoting.
    """
    total = len(accessions)
    steps: list[tuple[str, str]] = []
    for index, accession in enumerate(accessions, start=1):
        accession_dir = sra_dir / accession
        sra_file = accession_dir / f"{accession}.sra"
        fasta_file = accession_dir / f"{accession}.fasta"
        db_prefix = accession_dir / accession
        tag = f"[run {index}/{total}] {accession}"
        steps.append((
            f"{tag} - step 1/3: downloading .sra",
            f'"{prefetch}" --progress --max-size u -O "{sra_dir}" {accession}',
        ))
        steps.append((
            f"{tag} - step 2/3: converting to FASTA",
            f'"{fasterq_dump}" --fasta --split-spot --skip-technical --progress -f '
            f"--seq-defline '$ac.$si.$ri' "
            f'-O "{accession_dir}" -t "{accession_dir}" "{sra_file}"',
        ))
        steps.append((
            f"{tag} - step 3/3: building BLAST database",
            f'"{makeblastdb}" -in "{fasta_file}" -dbtype nucl -parse_seqids '
            f'-title {accession} -out "{db_prefix}"',
        ))
    return steps


def _windows_fetch_script(steps: list[tuple[str, str]]) -> str:
    """PowerShell body: hold the machine awake, run each step, release + pause.

    On Windows, *locking* the screen never stops a running process; only *sleep*
    does. So the fetch doesn't need `screen` (a Unix terminal-detach tool that
    can't touch power state anyway) -- it needs Windows told not to idle-sleep
    while it runs. SetThreadExecutionState(ES_CONTINUOUS|ES_SYSTEM_REQUIRED) does
    exactly that, scoped to this window and released in `finally`, so an
    unattended locked machine keeps downloading without changing the user's
    global power plan or needing admin. We *also* request away mode best-effort so
    a *manual* Sleep keeps the CPU running (screen off) where the machine supports
    it; a real S3 suspend powers the CPU/NIC off, so no code can download through
    it -- preventing the suspend is the only option. Each native step prints its progress
    header, then is guarded on $LASTEXITCODE so a failure stops the chain instead
    of building on a bad file. The guard runs *after* the native command, so
    $LASTEXITCODE is always set by then (the Write-Host lines never touch it).
    0x80000000/0x80000001 are cast from hex strings because PowerShell parses the
    bare literal 0x80000000 as a negative Int32 that won't convert to uint.
    """
    guard = "  if ($LASTEXITCODE -ne 0) { throw 'A step failed.' }\n"
    body = "".join(
        f"  Write-Host ''\n  Write-Host '=== {header} ==='\n  & {command}\n{guard}"
        for header, command in steps
    )
    return (
        "$p = Add-Type -PassThru -Name CoblastPwr -Namespace Win32 -MemberDefinition "
        "'[DllImport(\"kernel32.dll\")] public static extern uint SetThreadExecutionState(uint f);'\n"
        "[void]$p::SetThreadExecutionState([uint32]'0x80000001')  # ES_CONTINUOUS|ES_SYSTEM_REQUIRED: block idle sleep\n"
        # ponytail: best-effort away mode -- if the machine supports it and the power
        # plan enables it (desktops mostly; off by default; unsupported on most laptops
        # and Modern Standby), a manual Sleep keeps the CPU running with the screen off
        # so the download survives. Set AFTER the base lock: if away mode is unavailable
        # this call is a no-op and the idle-sleep block above still stands (no regression).
        # Ceiling: a real S3 suspend can't be worked around; that needs the lid/sleep-
        # button power setting changed, which we deliberately don't touch.
        "[void]$p::SetThreadExecutionState([uint32]'0x80000041')  # + ES_AWAYMODE_REQUIRED (best-effort)\n"
        "try {\n"
        + body
        + "  Write-Host ''\n"
        "  Write-Host '=== All runs downloaded, converted, and indexed. Refresh the SRA"
        " Workbench, then register the new BLASTDBs. ==='\n"
        "} catch {\n"
        "  Write-Host ''\n"
        "  Write-Host '*** A step above FAILED. Scroll up to read the error, fix it, then"
        " re-run. Nothing was registered. ***'\n"
        "} finally {\n"
        "  [void]$p::SetThreadExecutionState([uint32]'0x80000000')  # release the wake lock\n"
        "  Read-Host 'Press Enter to close'\n"
        "}\n"
    )


def _run_in_new_terminal(steps: list[tuple[str, str]]) -> None:
    """Write the (header, command) steps to a temp script and run it in a terminal.

    A temp script sidesteps cmd.exe/bash inline-quoting differences and lets one
    window run the whole prefetch+convert chain for every accession. Each step's
    progress header is echoed before its command so the user can watch position.
    """
    if os.name == "nt":
        fd, path = tempfile.mkstemp(prefix="coblast_sra_fetch_", suffix=".ps1")
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(_windows_fetch_script(steps))
        subprocess.Popen(
            ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", path],
            creationflags=subprocess.CREATE_NEW_CONSOLE,
        )
        return
    fd, path = tempfile.mkstemp(prefix="coblast_sra_fetch_", suffix=".sh")
    body = "".join(
        f"echo\necho '=== {header} ==='\n{command}\n" for header, command in steps
    )
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(
            "#!/bin/sh\nset -e\n"
            + body
            + "echo\necho '=== Done. Refresh the SRA Workbench and register the new BLASTDBs. ==='\n"
        )
    os.chmod(path, 0o755)
    if sys.platform == "darwin":
        run = json.dumps(f"sh {shlex.quote(path)}")
        subprocess.Popen(["osascript", "-e", f"tell application \"Terminal\" to do script {run}"])
    else:
        # ponytail: Linux terminal emulators vary too much to guess; run detached
        # in the background. Swap in `x-terminal-emulator -e sh path` for a window.
        subprocess.Popen(["/bin/sh", path], start_new_session=True)


def spawn_sra_fetch(accessions: list[str]) -> Path:
    """Fetch + convert the given runs in a new terminal; return the download dir.

    Runs in its own window so the web request returns immediately; the finished
    .sra and FASTA files appear in the SRA Projects table on the next page load.
    """
    sra_dir = sra_download_dir()
    steps = build_fetch_script_lines(
        accessions,
        sra_dir,
        sra_tool_exe("prefetch"),
        sra_tool_exe("fasterq-dump"),
        blast_exe("makeblastdb"),
    )
    _run_in_new_terminal(steps)
    return sra_dir


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
        "--split-spot",  # one record per biological read; no chimeric mate concat
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
    # Name from the database's own BLAST title (falling back to the prefix stem)
    # rather than a forced "SRA <acc> reads" wrapper. ponytail: one extra info
    # read; registration is a deliberate, infrequent action, not a hot path.
    info = verify_database_prefix(db_prefix_path)
    display_name = str(info.get("database_title") or "").strip() or Path(str(db_prefix_path)).name
    return register_existing_database(
        display_name=display_name,
        db_type="nucl",
        db_prefix_path=db_prefix_path,
        source_fasta_path=source,
        description=f"Local nucleotide BLAST database prepared from {accession}.",
        category="sra",
        notes="Registered from the SRA workbench.",
    )


if __name__ == "__main__":
    # ponytail: smallest check that fails if the accession guard breaks.
    assert parse_run_accessions("srr1, ERR2 ; DRR3\nSRR1") == ["SRR1", "ERR2", "DRR3"]
    assert parse_run_accessions(" srr9 ") == ["SRR9"]
    for bad in ["", "SRP123", "PRJNA1", "SRX9", "hello"]:
        try:
            parse_run_accessions(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected rejection for {bad!r}")
    print("parse_run_accessions OK")

    prefetch, fasterq_dump = Path("/tools/prefetch"), Path("/tools/fasterq-dump")
    makeblastdb = Path("/tools/makeblastdb")
    steps = build_fetch_script_lines(["SRR1", "SRR2"], Path("/data/sra"), prefetch, fasterq_dump, makeblastdb)
    headers = [header for header, _ in steps]
    cmds = [cmd for _, cmd in steps]
    assert len(steps) == 6  # prefetch + convert + index per accession
    assert cmds[0] == f'"{prefetch}" --progress --max-size u -O "{Path("/data/sra")}" SRR1'  # bar + no cap
    assert headers[0] == "[run 1/2] SRR1 - step 1/3: downloading .sra"  # position the user reads
    assert headers[3] == "[run 2/2] SRR2 - step 1/3: downloading .sra"  # second run counted
    assert str(fasterq_dump) in cmds[1] and "--progress" in cmds[1]  # native convert progress bar
    assert "--split-spot" in cmds[1]  # every read its own record, no chimeric mates
    assert "--seq-defline '$ac.$si.$ri'" in cmds[1]  # unique <acc>.<spot>.<read> ids, no mate collision
    assert str(makeblastdb) in cmds[2] and "-parse_seqids" in cmds[2]  # blast-ready, id index for eToL
    print("build_fetch_script_lines OK")

    win = _windows_fetch_script([("run 1", '"pf" -O "d" A'), ("run 2", '"fq" -x "d" A'),
                                 ("run 3", '"mb" -in "d" -out "p"')])
    assert "SetThreadExecutionState([uint32]'0x80000001')" in win  # holds the wake lock
    assert "SetThreadExecutionState([uint32]'0x80000041')" in win  # best-effort away mode
    assert win.count("SetThreadExecutionState([uint32]'0x80000000')") == 1  # releases it once
    assert win.count("if ($LASTEXITCODE -ne 0)") == 3  # every step guarded
    assert win.count("Write-Host '=== run ") == 3  # every step prints its progress header
    print("_windows_fetch_script OK")
