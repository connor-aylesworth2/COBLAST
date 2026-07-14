"""Command-line launcher for the COBLAST+ local interface.

This script prepares the runtime environment, checks BLAST+, optionally builds
or reuses a virtual environment, runs a smoke test, and starts the Flask app.
When packaged with PyInstaller it can also run as a standalone bundle.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import webbrowser

from config import is_frozen, resource_root as bundle_root, tool_name


APP_NAME = "COBLAST"
MIN_PYTHON = (3, 11)
REQUIRED_BLAST_TOOLS = (
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "makeblastdb",
    "blastdbcmd",
)


class LauncherError(RuntimeError):
    """Expected setup/startup failures that should be shown cleanly to users."""

    pass


@lru_cache(maxsize=1)
def project_root() -> Path:
    """Find the repository/app root from env, executable location, or script path."""
    env_root = os.environ.get("COBLAST_PROJECT_ROOT")
    candidates: list[Path] = []

    if env_root:
        candidates.append(Path(env_root))

    if is_frozen():
        exe_dir = Path(sys.executable).resolve().parent
        candidates.extend([exe_dir, exe_dir.parent])
    else:
        candidates.append(Path(__file__).resolve().parent)

    for candidate in candidates:
        root = candidate.expanduser().resolve()
        if (root / "app.py").exists() and (root / "requirements.txt").exists():
            return root

    # The searched list is included because most support issues are path-related.
    searched = "\n".join(f"  - {candidate}" for candidate in candidates)
    raise LauncherError(
        "Could not locate the COBLAST project folder. Run this launcher from the "
        "repository root, place run_COBLAST.exe in the repository root or dist "
        "folder, or set COBLAST_PROJECT_ROOT.\n\n"
        f"Searched:\n{searched}"
    )


def has_bundled_app() -> bool:
    """Detect whether the frozen executable contains the full Flask app bundle."""
    root = bundle_root()
    return (
        (root / "templates" / "index.html").exists()
        and (root / "static" / "styles.css").exists()
        and (root / "blast" / "bin" / tool_name("blastn")).exists()
    )


def standalone_data_dir() -> Path:
    """Choose the mutable data directory for a standalone executable.

    Delegates to config.runtime_data_dir() so the launcher and the bundled app
    agree on a single per-user location. This is called before COBLAST_DATA_DIR
    is set, so when frozen it resolves to the stable per-user data folder.
    """
    from config import runtime_data_dir

    return runtime_data_dir()


def _prompt_for_data_dir() -> Path | None:
    """Ask for the data location in a small window; return a validated dir or None.

    The window owns the input: an editable path field (prefilled with the default)
    plus a Browse button that opens File Explorer. The earlier version made the
    native folder dialog the *only* way in, so testers who lost that dialog — one
    reported it vanishing after a few seconds — had no way to give a path at all.
    Here the field always works, and Browse is just a convenience for populating it.

    Returns None when tkinter is unavailable/headless or the user cancels, so the
    caller keeps the default location rather than failing to launch.
    """
    try:
        import tkinter as tk
    except Exception:
        print(
            f"[{APP_NAME}] Folder picker unavailable; using the default data "
            "location. Change it later on the Settings page, or relaunch with "
            "--data-dir <path>.",
            flush=True,
        )
        return None

    from config import runtime_data_dir, validate_data_dir
    from folder_picker import ask_directory

    default_dir = runtime_data_dir()  # usually %LOCALAPPDATA%\COBLAST_data
    try:
        default_dir.mkdir(parents=True, exist_ok=True)  # so Browse opens here
    except OSError:
        pass

    chosen: Path | None = None
    root = None
    try:
        root = tk.Tk()
        root.title("COBLAST+ data location")
        root.attributes("-topmost", True)
        frame = tk.Frame(root, padx=16, pady=12)
        frame.pack(fill="both", expand=True)

        tk.Label(
            frame,
            justify="left",
            wraplength=560,
            text=(
                "COBLAST+ keeps its databases, SRA downloads, and results in one "
                "folder.\n\nUse the suggested folder, or point it at a drive with "
                "more free space. The path must have no spaces."
            ),
        ).pack(anchor="w")

        row = tk.Frame(frame)
        row.pack(fill="x", pady=(10, 4))
        entry = tk.Entry(row, width=60)
        entry.insert(0, str(default_dir))
        entry.pack(side="left", fill="x", expand=True)
        status = tk.Label(frame, fg="#b00020", justify="left", wraplength=560)

        def browse() -> None:
            try:
                picked = ask_directory(
                    "Choose where COBLAST+ stores its data "
                    "(databases, SRA downloads, results)",
                    entry.get().strip() or str(default_dir),
                    parent=root,
                )
            except RuntimeError as exc:
                status.config(text=f"{exc} Type the folder path instead.")
                status.pack(anchor="w")
                return
            if picked:
                entry.delete(0, tk.END)
                entry.insert(0, picked)

        def use_it() -> None:
            nonlocal chosen
            try:
                # Validates spaces/filesystem here (a makeblastdb probe), so a bad
                # drive fails now instead of mid-run. Briefly freezes the window.
                chosen = validate_data_dir(entry.get())
            except ValueError as exc:
                status.config(text=str(exc))
                status.pack(anchor="w")
                return
            root.destroy()

        tk.Button(row, text="Browse...", command=browse).pack(side="left", padx=(8, 0))
        buttons = tk.Frame(frame)
        buttons.pack(anchor="e", pady=(10, 0))
        tk.Button(buttons, text="Cancel", command=root.destroy).pack(side="right")
        tk.Button(buttons, text="Use this folder", command=use_it, default="active").pack(
            side="right", padx=(0, 8)
        )
        entry.bind("<Return>", lambda _event: use_it())

        root.mainloop()
    except Exception:
        # A swallowed traceback is why the tester report had nothing to go on.
        traceback.print_exc()
        return None
    finally:
        if root is not None:
            try:
                root.destroy()
            except Exception:
                pass  # already destroyed by Use/Cancel
    return chosen


def resolve_data_dir(args: argparse.Namespace) -> Path:
    """Pick the data dir from --data-dir, the saved pointer, or a first-run picker."""
    from config import load_saved_data_dir, runtime_data_dir, save_data_dir

    if args.data_dir:
        try:
            return save_data_dir(args.data_dir)
        except ValueError as exc:
            raise LauncherError(str(exc)) from exc

    want_picker = args.pick_data_dir or (
        is_frozen()
        and load_saved_data_dir() is None
        and not os.environ.get("COBLAST_DATA_DIR")
    )
    if want_picker:
        picked = _prompt_for_data_dir()
        if picked is not None:
            return save_data_dir(picked)

    return runtime_data_dir()  # env -> saved pointer -> frozen default -> source instance


def setup_data_location(args: argparse.Namespace) -> None:
    """Pin the data dir and redirect temp scratch onto the same drive.

    Exports COBLAST_DATA_DIR so every module and subprocess agrees on one
    location, and points TMP/TEMP at ``<data_dir>/tmp`` so large fasterq-dump /
    CAP3 / BLAST scratch does not quietly fill C:. Setting os.environ in the
    parent covers both the frozen in-process app and the source-mode child
    (which copies os.environ).
    """
    data_dir = resolve_data_dir(args)
    # A saved pointer or env value skips validate_data_dir, so re-check here that
    # BLAST+ can actually build databases on this filesystem (exFAT/FAT USB sticks
    # pass a plain write test but fail every makeblastdb). Clear message beats the
    # raw NCBI "no write permissions" exception deep in the first run.
    from config import require_blast_capable_data_dir

    try:
        require_blast_capable_data_dir(data_dir)
    except ValueError as exc:
        raise LauncherError(str(exc)) from exc
    data_dir.mkdir(parents=True, exist_ok=True)
    os.environ["COBLAST_DATA_DIR"] = str(data_dir)

    temp_dir = data_dir / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    os.environ["TMP"] = os.environ["TEMP"] = os.environ["TMPDIR"] = str(temp_dir)
    tempfile.tempdir = None  # drop any cached temp dir so the new TMP/TEMP takes effect


def step(message: str) -> None:
    """Print a visible progress heading in the terminal."""
    print(f"\n[{APP_NAME}] {message}", flush=True)


def display_command(command: list[str]) -> str:
    """Format a command for humans without changing the command that runs."""
    if os.name == "nt":
        return subprocess.list2cmdline(command)
    return " ".join(command)


def run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a setup command and turn non-zero exits into LauncherError."""
    printable = display_command(command)
    print(f"  > {printable}", flush=True)
    completed = subprocess.run(
        command,
        cwd=str(cwd or project_root()),
        env=env,
        text=True,
        check=False,
    )
    if check and completed.returncode != 0:
        raise LauncherError(f"Command failed with exit code {completed.returncode}: {printable}")
    return completed


def port_is_available(port: int) -> bool:
    """Check whether localhost can bind to the requested port."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        try:
            sock.bind(("127.0.0.1", port))
        except OSError:
            return False
    return True


def choose_available_port(requested_port: int) -> int:
    """Use the requested port or the next nearby free localhost port."""
    if port_is_available(requested_port):
        return requested_port

    for port in range(requested_port + 1, min(requested_port + 100, 65535) + 1):
        if port_is_available(port):
            print(
                f"  Port {requested_port} is already in use; using 127.0.0.1:{port} instead.",
                flush=True,
            )
            return port

    raise LauncherError(
        f"Could not find an available localhost port near {requested_port}. "
        "Close any existing COBLAST/Flask windows and try again."
    )


def open_browser_url(url: str) -> None:
    """Try to open the local app URL, falling back to terminal instructions."""
    try:
        opened = webbrowser.open(url)
    except Exception as exc:
        print(
            f"  Could not open the browser automatically: {exc}\n"
            f"  Open this address manually instead: {url}",
            flush=True,
        )
        return

    if not opened:
        print(
            "  Could not open the browser automatically.\n"
            f"  Open this address manually instead: {url}",
            flush=True,
        )


def wait_for_port(port: int, timeout: float = 30.0) -> bool:
    """Block until 127.0.0.1:port accepts a connection, or timeout elapses.

    A fixed sleep raced Flask's bind on slow/AV-scanned first runs, so the
    browser opened before the server was listening and testers saw a dead
    'port not accessible' tab. Polling the real socket removes the guess.
    """
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.5):
                return True
        except OSError:
            time.sleep(0.1)  # refused/unreachable yet; avoid a busy loop
    return False


def open_browser_later(url: str, port: int) -> None:
    """Open the browser once Flask is actually listening on port."""
    def opener() -> None:
        if not wait_for_port(port):
            print(
                "  The local server did not start listening in time.\n"
                f"  Open this address manually once it does: {url}",
                flush=True,
            )
            return
        open_browser_url(url)

    thread = threading.Thread(target=opener, daemon=True)
    thread.start()


def require_supported_python() -> None:
    """Ensure the current launcher process is running on supported Python."""
    if sys.version_info < MIN_PYTHON:
        version = ".".join(str(part) for part in MIN_PYTHON)
        current = ".".join(str(part) for part in sys.version_info[:3])
        raise LauncherError(
            f"{APP_NAME} needs Python {version} or newer. "
            f"This script is running with Python {current} at {sys.executable}."
        )


def python_probe(command: list[str]) -> tuple[Path, tuple[int, int, int]] | None:
    """Return executable path/version when a Python command is usable."""
    try:
        completed = subprocess.run(
            command
            + [
                "-c",
                (
                    "import sys; "
                    "print(sys.executable); "
                    "print(f'{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}')"
                ),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    if completed.returncode != 0:
        return None

    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        return None

    try:
        version = tuple(int(part) for part in lines[1].split(".")[:3])
    except ValueError:
        return None

    if len(version) != 3 or version < (*MIN_PYTHON, 0):
        return None

    return Path(lines[0]), version


def candidate_python_commands(cli_python: str | None) -> list[list[str]]:
    """Build an ordered, de-duplicated list of Python commands to try."""
    candidates: list[list[str]] = []

    if cli_python:
        candidates.append([cli_python])

    env_python = os.environ.get("COBLAST_PYTHON") or os.environ.get("PYTHON")
    if env_python:
        candidates.append([env_python])

    if not is_frozen():
        candidates.append([sys.executable])

    if os.name == "nt":
        # The Windows py launcher is often the most reliable way to select 3.11+.
        candidates.extend(
            [
                ["py", "-3.13"],
                ["py", "-3.12"],
                ["py", "-3.11"],
                ["py", "-3"],
            ]
        )

    for name in ("python", "python3"):
        path = shutil.which(name)
        if path:
            candidates.append([path])

    unique: list[list[str]] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = "\0".join(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def find_host_python(cli_python: str | None) -> list[str]:
    """Find a host Python capable of creating the virtual environment."""
    searched: list[str] = []
    for command in candidate_python_commands(cli_python):
        searched.append(display_command(command))
        if python_probe(command) is not None:
            return command

    version = ".".join(str(part) for part in MIN_PYTHON)
    searched_text = "\n".join(f"  - {command}" for command in searched)
    raise LauncherError(
        f"Could not find a usable Python {version}+ installation to create .venv.\n"
        "Install Python, or rerun with:\n"
        "  run_COBLAST.exe --python \"C:\\Path\\To\\python.exe\"\n\n"
        f"Tried:\n{searched_text}"
    )


def bundled_blast_bin() -> Path:
    """Return the BLAST+ bin path inside a standalone bundle."""
    return bundle_root() / "blast" / "bin"


def candidate_blast_bins(cli_blast_bin: str | None) -> list[Path]:
    """Build an ordered, de-duplicated list of BLAST+ bin folders to try."""
    root = project_root()
    candidates: list[Path] = []

    if cli_blast_bin:
        candidates.append(Path(cli_blast_bin))

    env_blast_bin = os.environ.get("BLAST_BIN")
    if env_blast_bin:
        candidates.append(Path(env_blast_bin))

    candidates.extend(
        [
            root / "ncbi-blast-2.17.0+" / "bin",
            root.parent / "ncbi-blast-2.17.0+" / "bin",
            Path(r"C:\Program Files\NCBI\blast-2.17.0+\bin"),
            Path(r"C:\Program Files\NCBI\blast+\bin"),
        ]
    )

    path_blastn = shutil.which(tool_name("blastn"))
    if path_blastn:
        candidates.append(Path(path_blastn).parent)

    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        try:
            normalized = str(candidate.expanduser().resolve())
        except OSError:
            normalized = str(candidate.expanduser())
        if normalized not in seen:
            seen.add(normalized)
            unique.append(Path(normalized))
    return unique


def missing_blast_tools(blast_bin: Path) -> list[str]:
    """List required BLAST+ executables missing from a candidate bin folder."""
    return [
        tool_name(tool)
        for tool in REQUIRED_BLAST_TOOLS
        if not (blast_bin / tool_name(tool)).exists()
    ]


def find_blast_bin(cli_blast_bin: str | None) -> Path:
    """Find a BLAST+ bin folder containing every required executable.

    If none is found in a known location, fetch a portable BLAST+ into the
    per-user data dir so a fresh machine can launch without a manual install.
    """
    for candidate in candidate_blast_bins(cli_blast_bin):
        if candidate.exists() and not missing_blast_tools(candidate):
            return candidate

    install_error: Exception | None = None
    try:
        from config import blast_bin_dir, ensure_tool_bin

        step("BLAST+ not found locally; fetching a portable copy")
        installed = ensure_tool_bin("blast", blast_bin_dir)
        if installed is not None and not missing_blast_tools(installed):
            print(f"  Installed BLAST+ into {installed}")
            return installed
    except Exception as exc:  # fall through to the manual-install message below
        install_error = exc

    searched = "\n".join(f"  - {candidate}" for candidate in candidate_blast_bins(cli_blast_bin))
    hint = f"\nAuto-install failed: {install_error}\n" if install_error else "\n"
    raise LauncherError(
        "Could not find or install a complete BLAST+ bin directory."
        + hint
        + "Install/extract NCBI BLAST+ and rerun with:\n"
        "  python run_COBLAST.py --blast-bin \"C:\\Tools\\ncbi-blast-2.17.0+\\bin\"\n\n"
        f"Searched:\n{searched}"
    )


def verify_blast(blast_bin: Path, env: dict[str, str]) -> None:
    """Confirm BLAST+ can run before the web interface starts."""
    step("Checking BLAST+ executables")
    missing = missing_blast_tools(blast_bin)
    if missing:
        raise LauncherError(
            f"BLAST+ bin directory is missing required tool(s): {', '.join(missing)}"
        )

    blastn_path = blast_bin / tool_name("blastn")
    try:
        # Running blastn -version catches blocked/quarantined executables early.
        completed = subprocess.run(
            [str(blastn_path), "-version"],
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise LauncherError(
            "Windows could not run the BLAST+ executable bundled with COBLAST.\n"
            f"Tried: {blastn_path}\n\n"
            f"Original error: {exc}\n\n"
            "Common causes are antivirus quarantine, Windows SmartScreen, AppLocker, "
            "running from inside a ZIP file, or a protected/network-synced folder. "
            "Move the release folder to a normal local folder such as C:\\COBLAST, "
            "extract it fully, right-click COBLAST.exe > Properties > Unblock if "
            "that checkbox is present, then run again."
        ) from exc

    output = (completed.stdout or completed.stderr).strip()
    if completed.returncode != 0:
        raise LauncherError(f"Could not run blastn -version:\n{output}")
    print(f"  BLAST_BIN={blast_bin}")
    print(f"  {output.splitlines()[0] if output else 'blastn version detected'}")


def venv_dir() -> Path:
    """Location of the source-checkout virtual environment."""
    return project_root() / ".venv"


def venv_python() -> Path:
    """Return the platform-specific Python executable inside .venv."""
    if os.name == "nt":
        return venv_dir() / "Scripts" / "python.exe"
    return venv_dir() / "bin" / "python"


def python_is_usable(python_path: Path) -> bool:
    """Check whether an existing virtualenv Python still works."""
    return python_probe([str(python_path)]) is not None


def ensure_virtualenv(host_python: list[str] | None) -> Path:
    """Create or repair the local virtual environment, then return its Python."""
    step("Creating or reusing the Python virtual environment")
    python_path = venv_python()

    if python_path.exists() and not python_is_usable(python_path):
        # Virtualenvs can break when moved between machines or Python versions.
        print("  Existing .venv looks stale or broken; rebuilding it.")
        shutil.rmtree(venv_dir())

    if not python_path.exists():
        if host_python is None:
            raise LauncherError("A Python installation is required to create .venv.")
        run_command(host_python + ["-m", "venv", str(venv_dir())])

    if not python_path.exists():
        raise LauncherError(f"Virtual environment Python was not created at {python_path}")

    if not python_is_usable(python_path):
        raise LauncherError(f"Virtual environment Python is not usable at {python_path}")

    print(f"  Using {python_path}")
    return python_path


def install_requirements(python_path: Path, env: dict[str, str], skip_install: bool) -> None:
    """Install/update dependencies unless the caller opted out."""
    if skip_install:
        step("Skipping dependency installation")
        return

    requirements = project_root() / "requirements.txt"
    if not requirements.exists():
        raise LauncherError(f"Missing requirements file: {requirements}")

    step("Installing Python dependencies")
    run_command([str(python_path), "-m", "pip", "install", "--upgrade", "pip"], env=env)
    run_command([str(python_path), "-m", "pip", "install", "-r", str(requirements)], env=env)


def run_smoke_test(python_path: Path, env: dict[str, str], skip_smoke: bool) -> None:
    """Run the backend smoke test in source-checkout mode."""
    if skip_smoke:
        step("Skipping backend smoke test")
        return

    smoke_test = project_root() / "smoke_test.py"
    if not smoke_test.exists():
        raise LauncherError(f"Missing smoke test: {smoke_test}")

    step("Running backend smoke test")
    run_command([str(python_path), str(smoke_test)], env=env)


def run_standalone_smoke_test(skip_smoke: bool) -> None:
    """Run the backend smoke test inside a standalone bundle."""
    if skip_smoke:
        step("Skipping backend smoke test")
        return

    step("Running bundled backend smoke test")
    from smoke_test import main as smoke_main

    smoke_main()


def start_app(
    python_path: Path,
    env: dict[str, str],
    *,
    port: int,
    open_browser: bool,
    check_only: bool,
) -> None:
    """Start Flask from source mode using the prepared virtual environment."""
    if check_only:
        step("Check-only mode complete")
        return

    app = project_root() / "app.py"
    if not app.exists():
        raise LauncherError(f"Missing Flask app: {app}")

    port = choose_available_port(port)
    url = f"http://127.0.0.1:{port}"
    # The Flask app reads this via config.flask_port().
    env["BLAST_FLASK_PORT"] = str(port)

    step(f"Starting the local interface at {url}")
    process = subprocess.Popen([str(python_path), str(app)], env=env, cwd=str(project_root()))
    try:
        if open_browser:
            open_browser_later(url, port)
        print(f"\nCOBLAST is running at {url}. Press Ctrl+C in this terminal to stop it.")
        process.wait()
    except KeyboardInterrupt:
        print("\nStopping COBLAST...")
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


def start_standalone_app(
    *,
    port: int,
    open_browser: bool,
    check_only: bool,
) -> None:
    """Start Flask directly from the packaged executable process."""
    if check_only:
        step("Check-only mode complete")
        return

    port = choose_available_port(port)
    url = f"http://127.0.0.1:{port}"
    # In standalone mode the app runs in this process, so os.environ is enough.
    os.environ["BLAST_FLASK_PORT"] = str(port)

    step(f"Starting the bundled local interface at {url}")
    from app import app as flask_app
    from config import FLASK_HOST, flask_port

    if open_browser:
        open_browser_later(url, port)

    print(f"\nCOBLAST is running at {url}. Press Ctrl+C in this terminal to stop it.")
    flask_app.run(host=FLASK_HOST, port=flask_port(), debug=False)


def run_self_check() -> int:
    """Run the packaged end-to-end self-check and exit (build gate + diagnostic).

    Exercises the read-recovery, human-filter, and CAP3 stages that the toy smoke
    test skips, so a frozen build that silently drops them fails here instead of
    in front of a user. Runs in both source and frozen mode; when frozen it uses
    the bundled BLAST+/CAP3 exactly as a normal launch would, and a throwaway
    data dir so it never touches the user's registry.
    """
    step("Running packaged self-check (read recovery + human filter + CAP3)")
    if is_frozen():
        os.environ.setdefault("BLAST_BIN", str(bundled_blast_bin()))
    os.environ.setdefault(
        "COBLAST_DATA_DIR",
        str(Path(tempfile.gettempdir()) / "coblast_self_check_data"),
    )
    from frozen_self_check import run as run_check

    return run_check()


def run_standalone(args: argparse.Namespace) -> int:
    """Set up environment variables and run the bundled app path."""
    data_dir = standalone_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)
    os.environ["COBLAST_DATA_DIR"] = str(data_dir)

    if args.blast_bin:
        blast_bin = Path(args.blast_bin).expanduser().resolve()
    else:
        blast_bin = bundled_blast_bin()
    os.environ["BLAST_BIN"] = str(blast_bin)

    print(f"{APP_NAME} standalone bundle")
    print(f"  Bundle: {bundle_root()}")
    print(f"  Data:   {data_dir}")

    verify_blast(blast_bin, os.environ.copy())
    run_standalone_smoke_test(args.skip_smoke)
    start_standalone_app(
        port=args.port,
        open_browser=not args.no_browser,
        check_only=args.check_only,
    )
    return 0


def parse_args() -> argparse.Namespace:
    """Define command-line flags for setup, checks, and startup behavior."""
    parser = argparse.ArgumentParser(
        description=(
            "Set up the COBLAST virtual environment, verify BLAST+, run the smoke "
            "test, and start the local Flask interface."
        )
    )
    parser.add_argument(
        "--blast-bin",
        help="Path to the NCBI BLAST+ bin directory. Overrides BLAST_BIN.",
    )
    parser.add_argument(
        "--python",
        help=(
            "Path to a Python 3.11+ executable used to create .venv. "
            "Useful when running run_COBLAST.exe on a machine with multiple Python versions."
        ),
    )
    parser.add_argument(
        "--port",
        type=int,
        default=5000,
        help="Localhost port for the Flask interface. Default: 5000.",
    )
    parser.add_argument(
        "--skip-install",
        action="store_true",
        help="Do not install or update Python dependencies.",
    )
    parser.add_argument(
        "--skip-smoke",
        action="store_true",
        help="Start the app without running smoke_test.py first.",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not automatically open the browser.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Verify BLAST+, create/update the venv, and run the smoke test without starting Flask.",
    )
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="Run the packaged end-to-end self-check (read recovery + human filter "
        "+ CAP3) and exit. Used as a post-build gate against the frozen .exe.",
    )
    parser.add_argument(
        "--data-dir",
        help="Folder where COBLAST+ stores its data (databases, SRA downloads, "
        "results). Saved for future launches. Must be a space-free path.",
    )
    parser.add_argument(
        "--pick-data-dir",
        action="store_true",
        help="Show a folder picker to choose the COBLAST+ data location, even if "
        "one is already saved.",
    )
    return parser.parse_args()


def main() -> int:
    """Launcher workflow used by both script and executable entry points."""
    args = parse_args()
    try:
        if args.port < 1 or args.port > 65535:
            raise LauncherError("Port must be between 1 and 65535.")

        require_supported_python()
        if args.self_check:
            return run_self_check()

        # Resolve the data location (flag/saved pointer/first-run picker) and
        # redirect temp scratch before any data-dir-derived module imports.
        setup_data_location(args)

        if is_frozen() and has_bundled_app():
            return run_standalone(args)

        # Source mode needs a host Python to create .venv unless a good .venv
        # already exists from an earlier launch.
        host_python = (
            None
            if venv_python().exists() and python_is_usable(venv_python())
            else find_host_python(args.python)
        )
        blast_bin = find_blast_bin(args.blast_bin)
        env = os.environ.copy()
        env["BLAST_BIN"] = str(blast_bin)

        verify_blast(blast_bin, env)
        python_path = ensure_virtualenv(host_python)
        install_requirements(python_path, env, args.skip_install)
        run_smoke_test(python_path, env, args.skip_smoke)
        start_app(
            python_path,
            env,
            port=args.port,
            open_browser=not args.no_browser,
            check_only=args.check_only,
        )
        return 0
    except LauncherError as exc:
        print(f"\n{APP_NAME} setup could not continue:\n{exc}", file=sys.stderr)
        print(
            "\nTroubleshooting hints:\n"
            "  - On Windows, try: py -3.11 run_COBLAST.py\n"
            "  - Confirm Python is 3.11 or newer: python --version\n"
            "  - For the .exe, pass --python if Windows finds the wrong Python.\n"
            "  - Confirm BLAST+ is installed and pass --blast-bin if needed.\n"
            "  - If dependency installation fails, check your internet connection and rerun.",
            file=sys.stderr,
        )
        return 1
    except OSError as exc:
        print(f"\n{APP_NAME} could not start because Windows denied access:\n{exc}", file=sys.stderr)
        print(
            "\nTroubleshooting hints:\n"
            "  - Extract the release folder fully before running COBLAST.exe.\n"
            "  - Move the release folder to a local folder such as C:\\COBLAST.\n"
            "  - Right-click COBLAST.exe > Properties > Unblock, if Windows shows that checkbox.\n"
            "  - Rerun from PowerShell with: .\\COBLAST.exe --check-only --skip-smoke --no-browser\n"
            "  - If your organization manages the computer, ask IT whether unsigned apps or "
            "temporary bundled executables are blocked.",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
