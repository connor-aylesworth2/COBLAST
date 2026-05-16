from pathlib import Path
import os


DEFAULT_BLAST_BIN = Path(r"C:\Program Files\NCBI\blast-2.17.0+\bin")
FLASK_HOST = "127.0.0.1"
DEFAULT_FLASK_PORT = 5000
REMOTE_BLAST_ENABLED = False
DISALLOWED_BLAST_OPTIONS = {"-remote"}


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


def flask_port() -> int:
    raw_port = os.environ.get("BLAST_FLASK_PORT", str(DEFAULT_FLASK_PORT))
    try:
        port = int(raw_port)
    except ValueError as exc:
        raise ValueError("BLAST_FLASK_PORT must be a whole number.") from exc
    if port < 1 or port > 65535:
        raise ValueError("BLAST_FLASK_PORT must be between 1 and 65535.")
    return port
