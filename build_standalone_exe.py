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

# Local modules bundled as hidden imports. PyInstaller's import graph misses
# modules imported lazily inside functions (e.g. run_COBLAST's late imports), so
# name them explicitly. One list feeds both the existence check and the
# --hidden-import flags so the two cannot drift. run_COBLAST (the entry script)
# and config (imported normally, so PyInstaller finds it) are handled separately.
MODULES = [
    "app",
    "blast_runner",
    "database_registry",
    "result_store",
    "database_size",
    "human_filter",
    "sra_workflow",
    "smoke_test",
    "frozen_self_check",
    "apoe_summary",
    "assembler",
    "contig_id",
    "design_matrix",
    "etol_summary",
    "etol_validation",
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


def cap3_filename() -> str:
    """Return the CAP3 executable filename for this OS."""
    return "cap3.exe" if os.name == "nt" else "cap3"


def default_cap3_bin() -> Path | None:
    """Choose the default CAP3 bin folder, if one is configured.

    Unlike BLAST+, CAP3 has no conventional install path, so this only honours
    the ``CAP3_BIN`` environment variable and otherwise returns None (build
    without CAP3).
    """
    env_cap3_bin = os.environ.get("CAP3_BIN")
    if env_cap3_bin:
        return Path(env_cap3_bin).expanduser().resolve()
    return None


def cap3_binary_files(cap3_bin: Path) -> list[Path]:
    """Return the CAP3 executable and any sibling DLLs it needs, if present.

    The Windows CAP3 build ships alongside runtime DLLs (e.g. cygwin1.dll), so
    every ``.dll`` next to the executable is bundled too. Returns an empty list
    when no CAP3 executable is found under ``cap3_bin``.
    """
    exe = cap3_bin / cap3_filename()
    if not exe.exists():
        return []
    return [exe, *sorted(cap3_bin.glob("*.dll"))]


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


def build_command(blast_bin: Path, cap3_bin: Path | None, name: str) -> list[str]:
    """Assemble the PyInstaller command for the standalone executable."""
    root = project_root()
    workpath = Path(tempfile.gettempdir()) / f"coblast_pyinstaller_{os.getpid()}"
    required_paths = [
        root / "run_COBLAST.py",
        root / "config.py",
        *[root / f"{module}.py" for module in MODULES],
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
        # UPX compression is a top antivirus false-positive trigger; skip it so
        # testers hit fewer quarantine/SmartScreen scares on an unsigned exe.
        "--noupx",
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

    for module in MODULES:
        command += ["--hidden-import", module]

    for filename in REQUIRED_BLAST_FILES + OPTIONAL_BLAST_FILES:
        path = blast_bin / filename
        if path.exists():
            command.extend(["--add-binary", add_data_arg(path, "blast/bin")])

    # CAP3 is optional: bundle it under cap3/bin when a directory is supplied so
    # config.cap3_bin_dir() finds it in the frozen app. Without it the exe still
    # runs but skips contig assembly/re-probing (the callers gate on it).
    if cap3_bin is not None:
        cap3_files = cap3_binary_files(cap3_bin)
        if not cap3_files:
            print(
                f"Warning: no CAP3 executable found in {cap3_bin}; building "
                "without CAP3 (contig assembly/re-probing will be skipped)."
            )
        for path in cap3_files:
            command.extend(["--add-binary", add_data_arg(path, "cap3/bin")])
    else:
        print(
            "No --cap3-bin/CAP3_BIN supplied; building without CAP3 (contig "
            "assembly/re-probing will be skipped in the packaged app)."
        )

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
    default_cap3 = default_cap3_bin()
    parser.add_argument(
        "--cap3-bin",
        default=str(default_cap3) if default_cap3 else None,
        help="Path to a directory containing the CAP3 executable to bundle. "
        "Optional (defaults to CAP3_BIN); without it the app runs but skips "
        "contig assembly/re-probing.",
    )
    parser.add_argument(
        "--name",
        default="COBLAST",
        help="Executable name to create under dist/. Default: COBLAST.",
    )
    parser.add_argument(
        "--skip-self-check",
        action="store_true",
        help="Do not run the packaged self-check after building. The self-check "
        "runs the freshly built .exe against a synthetic sample to prove read "
        "recovery, human filtering, and CAP3 assembly survived bundling.",
    )
    return parser.parse_args()


def main() -> int:
    """Run PyInstaller and print the built executable path when successful."""
    args = parse_args()
    blast_bin = Path(args.blast_bin).expanduser().resolve()
    print(f"Bundling BLAST+ from: {blast_bin}")

    cap3_bin = Path(args.cap3_bin).expanduser().resolve() if args.cap3_bin else None
    if cap3_bin is not None:
        print(f"Bundling CAP3 from: {cap3_bin}")

    if shutil.which("pyinstaller") is None:
        print("PyInstaller is not on PATH; using python -m PyInstaller.")

    command = build_command(blast_bin, cap3_bin, args.name)
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

    # Post-build gate: drive the freshly built .exe through the packaged
    # self-check. This is the only step that runs the *frozen* binary against a
    # known-answer sample, so a bundling regression (broken read recovery, a
    # missing/broken CAP3) fails the build here instead of a user's analysis.
    if not args.skip_self_check:
        if not exe_path.exists():
            print("Cannot run self-check: built executable not found.", file=sys.stderr)
            return 1
        print(f"\nRunning packaged self-check: {exe_path} --self-check")
        check = subprocess.run([str(exe_path), "--self-check"], cwd=project_root(), check=False)
        if check.returncode != 0:
            print(
                "\nPackaged self-check FAILED; the build is not trustworthy. "
                "Re-run with --skip-self-check only if you understand why.",
                file=sys.stderr,
            )
            return check.returncode
        print("Packaged self-check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
