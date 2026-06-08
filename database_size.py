"""Local BLAST database size helpers.

These functions report the on-disk footprint of a registered BLAST database
prefix and format byte counts for display. (A dynamic runtime estimator used to
live here but was removed: its predictions were too inaccurate to be useful to
clinicians, who in practice only vary which databases/patients are included.)
"""

from __future__ import annotations

from pathlib import Path


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
