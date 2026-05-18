from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import webbrowser


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
    pass


def project_root() -> Path:
    return Path(__file__).resolve().parent


def step(message: str) -> None:
    print(f"\n[{APP_NAME}] {message}", flush=True)


def tool_name(name: str) -> str:
    return f"{name}.exe" if os.name == "nt" else name


def run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    printable = " ".join(command)
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


def require_supported_python() -> None:
    if sys.version_info < MIN_PYTHON:
        version = ".".join(str(part) for part in MIN_PYTHON)
        current = ".".join(str(part) for part in sys.version_info[:3])
        raise LauncherError(
            f"{APP_NAME} needs Python {version} or newer. "
            f"This script is running with Python {current} at {sys.executable}."
        )


def candidate_blast_bins(cli_blast_bin: str | None) -> list[Path]:
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
    return [
        tool_name(tool)
        for tool in REQUIRED_BLAST_TOOLS
        if not (blast_bin / tool_name(tool)).exists()
    ]


def find_blast_bin(cli_blast_bin: str | None) -> Path:
    for candidate in candidate_blast_bins(cli_blast_bin):
        if candidate.exists() and not missing_blast_tools(candidate):
            return candidate

    searched = "\n".join(f"  - {candidate}" for candidate in candidate_blast_bins(cli_blast_bin))
    raise LauncherError(
        "Could not find a complete BLAST+ bin directory.\n"
        "Install/extract NCBI BLAST+ and rerun with:\n"
        "  python run_COBLAST.py --blast-bin \"C:\\Tools\\ncbi-blast-2.17.0+\\bin\"\n\n"
        f"Searched:\n{searched}"
    )


def verify_blast(blast_bin: Path, env: dict[str, str]) -> None:
    step("Checking BLAST+ executables")
    missing = missing_blast_tools(blast_bin)
    if missing:
        raise LauncherError(
            f"BLAST+ bin directory is missing required tool(s): {', '.join(missing)}"
        )

    completed = subprocess.run(
        [str(blast_bin / tool_name("blastn")), "-version"],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    output = (completed.stdout or completed.stderr).strip()
    if completed.returncode != 0:
        raise LauncherError(f"Could not run blastn -version:\n{output}")
    print(f"  BLAST_BIN={blast_bin}")
    print(f"  {output.splitlines()[0] if output else 'blastn version detected'}")


def venv_dir() -> Path:
    return project_root() / ".venv"


def venv_python() -> Path:
    if os.name == "nt":
        return venv_dir() / "Scripts" / "python.exe"
    return venv_dir() / "bin" / "python"


def python_is_usable(python_path: Path) -> bool:
    completed = subprocess.run(
        [str(python_path), "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.returncode == 0


def ensure_virtualenv() -> Path:
    step("Creating or reusing the Python virtual environment")
    python_path = venv_python()

    if python_path.exists() and not python_is_usable(python_path):
        print("  Existing .venv looks stale or broken; rebuilding it.")
        shutil.rmtree(venv_dir())

    if not python_path.exists():
        run_command([sys.executable, "-m", "venv", str(venv_dir())])

    if not python_path.exists():
        raise LauncherError(f"Virtual environment Python was not created at {python_path}")

    if not python_is_usable(python_path):
        raise LauncherError(f"Virtual environment Python is not usable at {python_path}")

    print(f"  Using {python_path}")
    return python_path


def install_requirements(python_path: Path, env: dict[str, str], skip_install: bool) -> None:
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
    if skip_smoke:
        step("Skipping backend smoke test")
        return

    smoke_test = project_root() / "smoke_test.py"
    if not smoke_test.exists():
        raise LauncherError(f"Missing smoke test: {smoke_test}")

    step("Running backend smoke test")
    run_command([str(python_path), str(smoke_test)], env=env)


def start_app(
    python_path: Path,
    env: dict[str, str],
    *,
    port: int,
    open_browser: bool,
    check_only: bool,
) -> None:
    if check_only:
        step("Check-only mode complete")
        return

    app = project_root() / "app.py"
    if not app.exists():
        raise LauncherError(f"Missing Flask app: {app}")

    url = f"http://127.0.0.1:{port}"
    env["BLAST_FLASK_PORT"] = str(port)

    step(f"Starting the local interface at {url}")
    process = subprocess.Popen([str(python_path), str(app)], env=env, cwd=str(project_root()))
    try:
        if open_browser:
            time.sleep(2)
            webbrowser.open(url)
        print("\nCOBLAST is running. Press Ctrl+C in this terminal to stop it.")
        process.wait()
    except KeyboardInterrupt:
        print("\nStopping COBLAST...")
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


def parse_args() -> argparse.Namespace:
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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.port < 1 or args.port > 65535:
            raise LauncherError("Port must be between 1 and 65535.")

        require_supported_python()
        blast_bin = find_blast_bin(args.blast_bin)
        env = os.environ.copy()
        env["BLAST_BIN"] = str(blast_bin)

        verify_blast(blast_bin, env)
        python_path = ensure_virtualenv()
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
            "  - Confirm BLAST+ is installed and pass --blast-bin if needed.\n"
            "  - If dependency installation fails, check your internet connection and rerun.",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
