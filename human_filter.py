"""Secondary human-genome read filter for the microbial eToL presets.

After the first eToL BLAST (microbial probes vs a patient database) is reduced to
the net of probe hits, some of the matched patient reads can still be human-
derived (for example, a read whose probe-matching core is shared with the host).
This is unavoidable because sequence similarity is non-transitive: a read can be
human yet still match a probe that itself has no human match. This module
performs the second-round human filtering described in Hu, Haas & Lathe 2022: it
pulls the matched patient reads back out, BLASTs them against a human genome
database, and reports which reads score a strong human alignment so the caller
can drop those hits from the final eToL results.

Reads are recovered by their ``sseqid`` (which equals the FASTA record id):
first via ``blastdbcmd`` (works when the patient DB was built with
``-parse_seqids``), otherwise by scanning the database's stored source FASTA.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Any

from config import blast_exe
from database_registry import blast_safe_path


# A read is removed when its best human-genome alignment scores above
# HUMAN_BITSCORE_THRESHOLD bits; the bitscore is the sole criterion. Hu, Haas &
# Lathe 2022 used a per-dataset cutoff adjusted to each library's mean read size
# (>160 MSBB, >126 Rockefeller, >100 Miami); >150 is the value they applied to
# brain (and liver/skin) datasets, so it is COBLAST's default for brain samples.
# Because the cutoff tracks read length, libraries with very different mean read
# lengths may warrant a different value. The E-value is left permissive so it
# never pre-filters an alignment, leaving bitscore as the only gate (no coverage
# or E-value filtering happens in this second, human-genome search).
DEFAULT_HUMAN_EVALUE = "1e9"
HUMAN_BITSCORE_THRESHOLD = 150.0
DEFAULT_HUMAN_TIMEOUT_SECONDS = 1800


def _unique_read_ids(hits: list[dict[str, str]]) -> list[str]:
    """Return the distinct patient-read ids (sseqid) referenced by the hits."""
    seen: dict[str, None] = {}
    for hit in hits:
        read_id = hit.get("sseqid", "")
        if read_id and read_id not in seen:
            seen[read_id] = None
    return list(seen)


def extract_reads_via_blastdbcmd(
    db_prefix_path: str, read_ids: list[str]
) -> dict[str, str] | None:
    """Pull reads from a BLAST DB by id; returns None if the DB lacks an id index."""
    if not read_ids:
        return {}
    with tempfile.TemporaryDirectory(prefix="human_filter_ids_") as tmpdir:
        ids_path = Path(tmpdir) / "ids.txt"
        ids_path.write_text("\n".join(read_ids) + "\n", encoding="utf-8")
        completed = subprocess.run(
            [
                str(blast_exe("blastdbcmd")),
                "-db",
                blast_safe_path(db_prefix_path),
                "-entry_batch",
                str(ids_path),
                "-outfmt",
                "%f",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    # No id index (built without -parse_seqids): signal the caller to fall back.
    if "no accession" in (completed.stderr or "").lower():
        return None
    if completed.returncode != 0 and not completed.stdout.strip():
        return None
    return _parse_fasta_text(completed.stdout)


def extract_reads_from_fasta(
    source_fasta_path: str, read_ids: list[str]
) -> dict[str, str]:
    """Recover reads by id from a source FASTA in a single streaming pass."""
    needed = set(read_ids)
    if not needed:
        return {}
    source = Path(source_fasta_path)
    if not source.exists():
        return {}

    reads: dict[str, str] = {}
    current_id: str | None = None
    keep = False
    parts: list[str] = []
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current_id is not None and keep:
                    reads[current_id] = "".join(parts)
                    if len(reads) == len(needed):
                        return reads
                header = line[1:].strip()
                current_id = header.split()[0] if header else None
                keep = current_id in needed
                parts = []
            elif keep:
                parts.append(line.strip())
    if current_id is not None and keep:
        reads[current_id] = "".join(parts)
    return reads


def _parse_fasta_text(text: str) -> dict[str, str]:
    reads: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current_id is not None:
                reads[current_id] = "".join(parts)
            header = line[1:].strip()
            current_id = header.split()[0] if header else None
            parts = []
        elif current_id is not None:
            parts.append(line.strip())
    if current_id is not None:
        reads[current_id] = "".join(parts)
    return reads


def extract_reads(
    db_prefix_path: str, source_fasta_path: str, read_ids: list[str]
) -> tuple[dict[str, str], str]:
    """Recover patient reads by id, returning (reads, method-used)."""
    db_reads = extract_reads_via_blastdbcmd(db_prefix_path, read_ids) or {}
    reads = dict(db_reads)
    methods = ["blastdbcmd"] if db_reads else []

    missing_ids = [read_id for read_id in read_ids if read_id not in reads]
    if source_fasta_path and missing_ids:
        fasta_reads = extract_reads_from_fasta(source_fasta_path, missing_ids)
        if fasta_reads:
            reads.update(fasta_reads)
            methods.append("source_fasta")

    return reads, "+".join(methods) if methods else "none"


def find_human_read_ids(
    reads: dict[str, str],
    human_db_prefix_path: str,
    *,
    evalue: str = DEFAULT_HUMAN_EVALUE,
    bitscore_threshold: float = HUMAN_BITSCORE_THRESHOLD,
    timeout_seconds: int = DEFAULT_HUMAN_TIMEOUT_SECONDS,
    num_threads: int | str | None = None,
) -> set[str]:
    """Return read ids whose best human-genome alignment scores above
    ``bitscore_threshold`` bits (Hu, Haas & Lathe 2022 brain cutoff = 150)."""
    if not reads:
        return set()
    with tempfile.TemporaryDirectory(prefix="human_filter_q_") as tmpdir:
        query_path = Path(tmpdir) / "reads.fasta"
        query_path.write_text(
            "".join(f">{read_id}\n{sequence}\n" for read_id, sequence in reads.items()),
            encoding="utf-8",
        )
        # -max_target_seqs 1 keeps the single best subject, whose top HSP carries
        # the read's maximum bitscore -- exactly the value the threshold gates on.
        command = [
            str(blast_exe("blastn")),
            "-task",
            "megablast",
            "-query",
            str(query_path),
            "-db",
            blast_safe_path(human_db_prefix_path),
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            "1",
            "-outfmt",
            "6 qseqid bitscore",
        ]
        # This search of tens of thousands of matched reads is the eToL pipeline's
        # heaviest post-BLAST step; run it across the job's cores rather than one.
        if num_threads:
            command += ["-num_threads", str(num_threads)]
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "Human-genome BLAST failed: "
            + ((completed.stderr or completed.stdout).strip() or "unknown error")
        )
    human_ids: set[str] = set()
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        try:
            bitscore = float(fields[1])
        except ValueError:
            continue
        if bitscore > bitscore_threshold:
            human_ids.add(fields[0].strip())
    return human_ids


def filter_human_hits(
    hits: list[dict[str, str]],
    *,
    db_prefix_path: str,
    source_fasta_path: str,
    human_db_prefix_path: str,
    evalue: str = DEFAULT_HUMAN_EVALUE,
    bitscore_threshold: float = HUMAN_BITSCORE_THRESHOLD,
    timeout_seconds: int = DEFAULT_HUMAN_TIMEOUT_SECONDS,
    num_threads: int | str | None = None,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Drop hits whose patient read scores a strong human-genome alignment.

    A read is treated as human when its best human alignment exceeds
    ``bitscore_threshold`` bits. Returns the kept hits plus a stats dict
    describing what happened. Reads that cannot be recovered are conservatively
    kept (never dropped on a guess).
    """
    read_ids = _unique_read_ids(hits)
    stats: dict[str, Any] = {
        "reads_total": len(read_ids),
        "reads_checked": 0,
        "reads_unresolved": len(read_ids),
        "human_reads": 0,
        "hits_removed": 0,
        "method": "none",
        "note": "",
    }
    if not read_ids:
        return hits, stats

    reads, method = extract_reads(db_prefix_path, source_fasta_path, read_ids)
    stats["method"] = method
    stats["reads_checked"] = len(reads)
    stats["reads_unresolved"] = len(read_ids) - len(reads)

    if not reads:
        stats["note"] = (
            "Could not recover patient reads for this database (no id-indexed "
            "BLAST DB and no readable source FASTA); human filter skipped."
        )
        return hits, stats

    human_ids = find_human_read_ids(
        reads,
        human_db_prefix_path,
        evalue=evalue,
        bitscore_threshold=bitscore_threshold,
        timeout_seconds=timeout_seconds,
        num_threads=num_threads,
    )
    stats["human_reads"] = len(human_ids)

    kept = [hit for hit in hits if hit.get("sseqid", "") not in human_ids]
    stats["hits_removed"] = len(hits) - len(kept)
    if stats["reads_unresolved"]:
        stats["note"] = (
            f"{stats['reads_unresolved']} matched read(s) could not be recovered "
            "and were kept unfiltered."
        )
    return kept, stats
