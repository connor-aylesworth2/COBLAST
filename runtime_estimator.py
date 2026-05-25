"""Rough local runtime estimates for BLAST jobs.

These estimates are intentionally conservative and are meant for planning pilot
runs, not for scheduling clinical production work. They become useful when
compared with observed runtimes from this same machine.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PROGRAM_MULTIPLIERS = {
    "blastn": 1.0,
    "blastp": 0.8,
    "blastx": 2.0,
    "tblastn": 2.4,
}
PRESET_MULTIPLIERS = {
    "fast": 0.6,
    "standard": 1.0,
    "sensitive": 1.8,
}
BLAST_DATABASE_SUFFIXES = {
    ".nal",
    ".ndb",
    ".nhr",
    ".nin",
    ".njs",
    ".not",
    ".nsq",
    ".ntf",
    ".nto",
    ".pal",
    ".pdb",
    ".phr",
    ".pin",
    ".pjs",
    ".pot",
    ".psq",
    ".ptf",
    ".pto",
}


@dataclass(frozen=True)
class RuntimeEstimate:
    seconds: float
    low_seconds: float
    high_seconds: float
    note: str


def database_prefix_files(db_prefix_path: str | Path) -> list[Path]:
    """Return files that appear to belong to a BLAST database prefix."""
    prefix = Path(db_prefix_path)
    parent = prefix.parent
    if not parent.exists():
        return []

    prefix_name = prefix.name
    return [
        path
        for path in parent.glob(f"{prefix_name}*")
        if path.is_file()
        and path.name.startswith(prefix_name)
        and path.suffix.lower() in BLAST_DATABASE_SUFFIXES
    ]


def database_storage_bytes(db_prefix_path: str | Path) -> int:
    """Sum the local files that share a BLAST database prefix."""
    total = 0
    for path in database_prefix_files(db_prefix_path):
        try:
            total += path.stat().st_size
        except OSError:
            continue
    return total


def format_bytes(size_bytes: int) -> str:
    """Format byte counts for compact UI display."""
    size = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{int(size)} B"
        size /= 1024
    return f"{size:.1f} TB"


def format_duration(seconds: float | None) -> str:
    """Format a duration estimate without pretending it is exact."""
    if seconds is None:
        return "unknown"
    if seconds < 60:
        return f"{seconds:.0f} sec"
    minutes = seconds / 60
    if minutes < 60:
        return f"{minutes:.1f} min"
    hours = minutes / 60
    if hours < 48:
        return f"{hours:.1f} hr"
    days = hours / 24
    return f"{days:.1f} days"


def estimate_blast_runtime_seconds(
    *,
    program: str,
    query_total_length: int,
    database_bytes: int,
    sensitivity_preset: str,
) -> RuntimeEstimate | None:
    """Estimate BLAST runtime from query size, database size, and program type."""
    if query_total_length <= 0 or database_bytes <= 0:
        return None

    database_gb = database_bytes / 1_000_000_000
    query_kb = max(query_total_length / 1_000, 0.05)
    program_multiplier = PROGRAM_MULTIPLIERS.get(program, 1.0)
    preset_multiplier = PRESET_MULTIPLIERS.get(sensitivity_preset, 1.0)

    # Local process startup and disk latency dominate tiny jobs; larger jobs
    # scale roughly with database size, query size, and translated-search cost.
    seconds = 2.0 + (database_gb * query_kb * program_multiplier * preset_multiplier * 90.0)
    low_seconds = max(1.0, seconds / 3)
    high_seconds = seconds * 3
    return RuntimeEstimate(
        seconds=seconds,
        low_seconds=low_seconds,
        high_seconds=high_seconds,
        note="rough estimate; calibrate with observed runtimes on this machine",
    )
