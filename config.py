from pathlib import Path
import os
import sys


DEFAULT_BLAST_BIN = Path(r"C:\Program Files\NCBI\blast-2.17.0+\bin")
FLASK_HOST = "127.0.0.1"
DEFAULT_FLASK_PORT = 5000
REMOTE_BLAST_ENABLED = False
DISALLOWED_BLAST_OPTIONS = {"-remote"}


def is_frozen() -> bool:
    return bool(getattr(sys, "frozen", False))


def resource_root() -> Path:
    if is_frozen():
        return Path(getattr(sys, "_MEIPASS")).resolve()
    return Path(__file__).resolve().parent


def resource_path(*parts: str) -> Path:
    return resource_root().joinpath(*parts)


def runtime_data_dir() -> Path:
    env_data_dir = os.environ.get("COBLAST_DATA_DIR")
    if env_data_dir:
        return Path(env_data_dir).expanduser().resolve()
    if is_frozen():
        return Path(sys.executable).resolve().parent / "COBLAST_data"
    return resource_root() / "instance"


def blast_bin_dir() -> Path:
    env_blast_bin = os.environ.get("BLAST_BIN")
    if env_blast_bin:
        return Path(env_blast_bin)

    bundled_blast_bin = resource_path("blast", "bin")
    if bundled_blast_bin.exists():
        return bundled_blast_bin

    return DEFAULT_BLAST_BIN


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
