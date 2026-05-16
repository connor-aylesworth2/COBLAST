from pathlib import Path
import os


DEFAULT_BLAST_BIN = Path(r"C:\Program Files\NCBI\blast-2.17.0+\bin")


def blast_bin_dir() -> Path:
    return Path(os.environ.get("BLAST_BIN", DEFAULT_BLAST_BIN))


def blast_exe(name: str) -> Path:
    suffix = ".exe" if os.name == "nt" and not name.endswith(".exe") else ""
    exe_path = blast_bin_dir() / f"{name}{suffix}"
    if not exe_path.exists():
        raise FileNotFoundError(
            f"Could not find {exe_path}. Set BLAST_BIN to your BLAST+ bin directory."
        )
    return exe_path
