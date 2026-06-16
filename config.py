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


def tool_name(name: str) -> str:
    """Add the Windows ``.exe`` suffix to a bare executable name when needed."""
    return f"{name}.exe" if os.name == "nt" and not name.endswith(".exe") else name


def blast_exe(name: str) -> Path:
    """Resolve one BLAST+ executable and fail with a useful setup message."""
    exe_path = blast_bin_dir() / tool_name(name)
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


# Cores held back from a default BLAST run so the machine stays responsive.
CPU_RESERVE_SMALL = 1  # machines with <= 4 logical cores
CPU_RESERVE_LARGE = 2  # machines with > 4 logical cores


def available_cpu_count() -> int:
    """Return the number of logical CPUs visible to this process (>= 1)."""
    return os.cpu_count() or 1


def default_thread_count() -> int:
    """Choose a sensible default for BLAST ``-num_threads``.

    Uses most of the machine but reserves a core or two so a clinician's
    desktop stays responsive during a search. The ``COBLAST_NUM_THREADS``
    environment variable and the per-job advanced field can override this.
    """
    total = available_cpu_count()
    reserve = CPU_RESERVE_SMALL if total <= 4 else CPU_RESERVE_LARGE
    return max(1, total - reserve)


COBLAST_BATCH_WORKERS_ENV = "COBLAST_BATCH_WORKERS"


def allocate_batch_resources(
    num_jobs: int, requested_workers: int | str | None = None
) -> tuple[int, int]:
    """Split the core budget into (concurrent workers, threads per job).

    Benchmarks show concurrency across patient databases scales far better than
    ``-num_threads`` within one search, so the budget is spent on workers first;
    leftover cores go to per-job threads only when there are fewer databases
    than the budget. Left on auto the product never exceeds the budget, so runs
    do not oversubscribe the CPU. Precedence for the worker count: an explicit
    ``requested_workers`` (the batch advanced field) > ``COBLAST_BATCH_WORKERS``
    > auto. A forced value that oversubscribes is the caller's responsibility.
    """
    budget = default_thread_count()
    jobs = max(1, int(num_jobs))
    workers = min(budget, jobs)

    requested = requested_workers
    if requested in (None, ""):
        requested = os.environ.get(COBLAST_BATCH_WORKERS_ENV)
    if requested not in (None, ""):
        try:
            workers = max(1, min(int(requested), jobs))
        except (TypeError, ValueError):
            pass

    threads_per_job = max(1, budget // workers)
    return workers, threads_per_job
