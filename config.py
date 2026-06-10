"""Shared configuration helpers for source and bundled COBLAST+ runs."""

from pathlib import Path
import os
import sys


# Default install location used when BLAST_BIN is not supplied.
DEFAULT_BLAST_BIN = Path(r"C:\Program Files\NCBI\blast-2.17.0+\bin")
FLASK_HOST = "127.0.0.1"
DEFAULT_FLASK_PORT = 5000
REMOTE_BLAST_ENABLED = False
DISALLOWED_BLAST_OPTIONS = {"-remote"}


def is_frozen() -> bool:
    """Return True when the app is running from a PyInstaller executable."""
    return bool(getattr(sys, "frozen", False))


def resource_root() -> Path:
    """Locate bundled/static resources in both source and frozen modes."""
    if is_frozen():
        return Path(getattr(sys, "_MEIPASS")).resolve()
    return Path(__file__).resolve().parent


def resource_path(*parts: str) -> Path:
    """Build a path inside the resource root."""
    return resource_root().joinpath(*parts)


# Folder name used for mutable per-user data in both source and frozen runs.
DATA_DIR_NAME = "COBLAST_data"


def user_data_base() -> Path:
    """Return the OS-specific base directory for per-user application data."""
    if os.name == "nt":
        local_app_data = os.environ.get("LOCALAPPDATA")
        if local_app_data:
            return Path(local_app_data)
        return Path.home() / "AppData" / "Local"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support"
    xdg_data_home = os.environ.get("XDG_DATA_HOME")
    if xdg_data_home:
        return Path(xdg_data_home)
    return Path.home() / ".local" / "share"


def runtime_data_dir() -> Path:
    """Choose where mutable app data should live."""
    env_data_dir = os.environ.get("COBLAST_DATA_DIR")
    if env_data_dir:
        return Path(env_data_dir).expanduser().resolve()
    if is_frozen():
        # A stable per-user location keeps data across version installs instead
        # of trapping it beside whichever .exe happened to create it.
        return (user_data_base() / DATA_DIR_NAME).resolve()
    return resource_root() / "instance"


def blast_bin_dir() -> Path:
    """Find the BLAST+ bin directory from env, bundle, or default install path."""
    env_blast_bin = os.environ.get("BLAST_BIN")
    if env_blast_bin:
        return Path(env_blast_bin)

    bundled_blast_bin = resource_path("blast", "bin")
    if bundled_blast_bin.exists():
        return bundled_blast_bin

    return DEFAULT_BLAST_BIN


def blast_exe(name: str) -> Path:
    """Resolve one BLAST+ executable and fail with a useful setup message."""
    suffix = ".exe" if os.name == "nt" and not name.endswith(".exe") else ""
    exe_path = blast_bin_dir() / f"{name}{suffix}"
    if not exe_path.exists():
        raise FileNotFoundError(
            f"Could not find {exe_path}. Set BLAST_BIN to your BLAST+ bin directory."
        )
    return exe_path


def flask_port() -> int:
    """Read and validate the Flask port from BLAST_FLASK_PORT."""
    raw_port = os.environ.get("BLAST_FLASK_PORT", str(DEFAULT_FLASK_PORT))
    try:
        port = int(raw_port)
    except ValueError as exc:
        raise ValueError("BLAST_FLASK_PORT must be a whole number.") from exc
    if port < 1 or port > 65535:
        raise ValueError("BLAST_FLASK_PORT must be between 1 and 65535.")
    return port
