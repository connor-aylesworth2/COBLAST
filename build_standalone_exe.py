"""Build a PyInstaller standalone COBLAST+ executable.

The generated executable bundles the Flask source files, sample data, and the
required BLAST+ binaries so the interface can start on machines without a local
checkout.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


REQUIRED_BLAST_FILES = [
    # These executables are necessary for every supported search/database path.
    "blastn.exe",
    "blastp.exe",
    "blastx.exe",
    "tblastn.exe",
    "makeblastdb.exe",
    "blastdbcmd.exe",
]
OPTIONAL_BLAST_FILES = [
    # Include supporting files when they exist, but do not fail older installs.
    "blastn.exe.manifest",
    "blastp.exe.manifest",
    "blastx.exe.manifest",
    "tblastn.exe.manifest",
    "makeblastdb.exe.manifest",
    "blastdbcmd.exe.manifest",
    "ncbi-vdb-md.dll",
    "nghttp2.dll",
]


def project_root() -> Path:
    """Return the repository root that contains this build script."""
    return Path(__file__).resolve().parent


def default_blast_bin() -> Path:
    """Choose the default BLAST+ bin folder used by the build command."""
    env_blast_bin = os.environ.get("BLAST_BIN")
    if env_blast_bin:
        return Path(env_blast_bin).expanduser().resolve()
    return (project_root().parent / "ncbi-blast-2.17.0+" / "bin").resolve()


def check_file(path: Path) -> None:
    """Fail early when a required source or binary file is missing."""
    if not path.exists():
        raise FileNotFoundError(path)


def pyinstaller_separator() -> str:
    """Return the source/destination separator used by PyInstaller on this OS."""
    return ";" if os.name == "nt" else ":"


def add_data_arg(source: Path, destination: str) -> str:
    """Format a PyInstaller --add-data/--add-binary argument."""
    return f"{source}{pyinstaller_separator()}{destination}"


def build_command(blast_bin: Path, name: str) -> list[str]:
    """Assemble the PyInstaller command for the standalone executable."""
    root = project_root()
    workpath = Path(tempfile.gettempdir()) / f"coblast_pyinstaller_{os.getpid()}"
    required_paths = [
        root / "run_COBLAST.py",
        root / "app.py",
        root / "blast_runner.py",
        root / "config.py",
        root / "database_registry.py",
        root / "result_store.py",
        root / "database_size.py",
        root / "human_filter.py",
        root / "sra_workflow.py",
        root / "smoke_test.py",
        root / "templates",
        root / "static",
        root / "sample_data",
        root / "data",
        root / "requirements.txt",
    ]
    for path in required_paths:
        check_file(path)

    for filename in REQUIRED_BLAST_FILES:
        check_file(blast_bin / filename)

    # A temp workpath keeps PyInstaller intermediate files out of the repo.
    command = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--clean",
        "--onefile",
        "--name",
        name,
        "--workpath",
        str(workpath),
        "--distpath",
        str(root / "dist"),
        "--collect-submodules",
        "Bio.SeqIO",
        "--exclude-module",
        "tkinter",
        "--exclude-module",
        "_tkinter",
        "--hidden-import",
        "app",
        "--hidden-import",
        "smoke_test",
        "--hidden-import",
        "blast_runner",
        "--hidden-import",
        "database_registry",
        "--hidden-import",
        "result_store",
        "--hidden-import",
        "database_size",
        "--hidden-import",
        "human_filter",
        "--hidden-import",
        "sra_workflow",
        "--add-data",
        add_data_arg(root / "templates", "templates"),
        "--add-data",
        add_data_arg(root / "static", "static"),
        "--add-data",
        add_data_arg(root / "sample_data", "sample_data"),
        "--add-data",
        add_data_arg(root / "data", "data"),
        "--add-data",
        add_data_arg(root / "requirements.txt", "."),
    ]

    for filename in REQUIRED_BLAST_FILES + OPTIONAL_BLAST_FILES:
        path = blast_bin / filename
        if path.exists():
            command.extend(["--add-binary", add_data_arg(path, "blast/bin")])

    command.append(str(root / "run_COBLAST.py"))
    return command


def parse_args() -> argparse.Namespace:
    """Read build-time options for BLAST+ location and executable name."""
    parser = argparse.ArgumentParser(
        description="Build a standalone COBLAST Windows executable with bundled BLAST+."
    )
    parser.add_argument(
        "--blast-bin",
        default=str(default_blast_bin()),
        help="Path to the NCBI BLAST+ bin directory to bundle.",
    )
    parser.add_argument(
        "--name",
        default="COBLAST",
        help="Executable name to create under dist/. Default: COBLAST.",
    )
    return parser.parse_args()


def main() -> int:
    """Run PyInstaller and print the built executable path when successful."""
    args = parse_args()
    blast_bin = Path(args.blast_bin).expanduser().resolve()
    print(f"Bundling BLAST+ from: {blast_bin}")

    if shutil.which("pyinstaller") is None:
        print("PyInstaller is not on PATH; using python -m PyInstaller.")

    command = build_command(blast_bin, args.name)
    print("Running:")
    print(" ".join(command))
    completed = subprocess.run(command, cwd=project_root(), check=False)
    if completed.returncode != 0:
        return completed.returncode

    exe_path = project_root() / "dist" / f"{args.name}.exe"
    if exe_path.exists():
        size_mb = exe_path.stat().st_size / (1024 * 1024)
        print(f"Built: {exe_path}")
        print(f"Size:  {size_mb:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
