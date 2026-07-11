"""Shared configuration helpers for source and bundled COBLAST+ runs."""

from pathlib import Path
import ctypes
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
import zipfile


# Default install location used when BLAST_BIN is not supplied.
DEFAULT_BLAST_BIN = Path(r"C:\Program Files\NCBI\blast-2.17.0+\bin")
FLASK_HOST = "127.0.0.1"
DEFAULT_FLASK_PORT = 5000
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


def _settings_path() -> Path:
    """Fixed-location pointer file recording the user's chosen data dir.

    Lives under the OS per-user base (on Windows: %LOCALAPPDATA%\\COBLAST), never
    inside the movable data dir it points at — that chicken-and-egg is why this
    tiny JSON exists instead of storing the choice in the SQLite registry (which
    itself lives in the data dir).
    """
    return user_data_base() / "COBLAST" / "settings.json"


def blast_incompatible_filesystem(path: str | Path) -> str | None:
    """Return the filesystem name if `path` is on one BLAST+ can't build DBs on.

    makeblastdb refuses to create databases on FAT/exFAT volumes (common on USB
    sticks and SD cards): its CreateDirectories() write-permission check needs the
    ACL model those filesystems lack, so it fails with "You do not have write
    permissions" even though the OS can write there fine. NTFS/ReFS are fine.
    Windows-only; returns None elsewhere (no drive volumes to inspect) and None
    when the volume can't be queried (e.g. an unmounted drive), so this never
    blocks on ambiguity — it only rejects a filesystem it positively identifies.
    """
    # ponytail: name-based FAT-family check (instant, no BLAST needed). If some
    # other filesystem ever shows the same makeblastdb failure, swap this for a
    # real one-shot makeblastdb probe of the directory.
    if os.name != "nt":
        return None
    drive = os.path.splitdrive(os.fspath(Path(path)))[0]
    if not drive:
        return None
    buf = ctypes.create_unicode_buffer(256)
    ok = ctypes.windll.kernel32.GetVolumeInformationW(
        ctypes.c_wchar_p(drive + "\\"), None, 0, None, None, None, buf, ctypes.sizeof(buf)
    )
    if not ok:
        return None
    return buf.value if buf.value.upper() in {"FAT", "FAT32", "EXFAT"} else None


def require_blast_capable_data_dir(path: str | Path) -> None:
    """Raise ValueError if BLAST+ cannot create databases on `path`'s filesystem.

    Fast filesystem-name check only (no BLAST+ process) so it is cheap to call on
    every launch. `makeblastdb_probe_error` is the authoritative version used when
    a new write-to location is chosen.
    """
    bad_fs = blast_incompatible_filesystem(path)
    if bad_fs:
        raise ValueError(
            f"The data folder {Path(path)} is on a {bad_fs} filesystem. BLAST+ cannot "
            "create databases on FAT/exFAT drives (common on USB sticks and SD cards). "
            "Choose a folder on an NTFS drive, or reformat the drive as NTFS."
        )


def _resolve_makeblastdb_no_install() -> Path | None:
    """Find makeblastdb from env/bundle/default install, WITHOUT auto-installing.

    Mirrors blast_bin_dir()'s precedence but never triggers a download, so using
    it to validate a directory can't kick off a multi-hundred-MB BLAST+ fetch.
    """
    candidates: list[Path] = []
    env_bin = os.environ.get("BLAST_BIN")
    if env_bin:
        candidates.append(Path(env_bin))
    candidates.append(resource_path("blast", "bin"))
    candidates.append(DEFAULT_BLAST_BIN)
    for directory in candidates:
        exe = directory / tool_name("makeblastdb")
        if exe.exists():
            return exe
    return None


def makeblastdb_probe_error(directory: str | Path) -> str | None:
    """Actually build a throwaway 1-sequence DB in `directory`; return an error or None.

    This is the authoritative "can BLAST+ write a database here" check that the
    fast filesystem-name check only approximates — it also catches network shares
    and other mounts that report an ordinary filesystem yet still reject
    makeblastdb. Returns None (skip) when makeblastdb is not already available, so
    validating a directory never forces a BLAST+ download.
    """
    makeblastdb = _resolve_makeblastdb_no_install()
    if makeblastdb is None:
        return None
    probe_dir = Path(directory) / ".coblast_blastdb_probe"
    try:
        probe_dir.mkdir(parents=True, exist_ok=True)
        fasta = probe_dir / "probe.fasta"
        fasta.write_text(">probe\nACGTACGTACGTACGT\n", encoding="utf-8")
        completed = subprocess.run(
            [str(makeblastdb), "-in", str(fasta), "-dbtype", "nucl", "-out", str(probe_dir / "probe")],
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            message = (completed.stderr or completed.stdout).strip()
            return message.splitlines()[-1].strip() if message else "makeblastdb failed."
        return None
    except OSError as exc:
        return str(exc)
    finally:
        shutil.rmtree(probe_dir, ignore_errors=True)


def _validate_writable_dir(path: str | Path) -> Path:
    """Resolve `path`, reject spaces, and confirm it is creatable/writable.

    Rejects paths with spaces because BLAST+ splits ``-out``/``-db`` on internal
    whitespace, so a spaced path silently yields broken databases. Shared by the
    data dir (write-to) and SRA reads dir (pull-from) validators.
    """
    raw = str(path).strip()
    if not raw:
        raise ValueError("Enter a folder path.")
    resolved = Path(raw).expanduser().resolve()
    if " " in str(resolved):
        raise ValueError(
            f"The folder path cannot contain spaces (got: {resolved}). BLAST+ cannot "
            r"build databases under a spaced path — choose a space-free folder such as D:\COBLAST_data."
        )
    try:
        resolved.mkdir(parents=True, exist_ok=True)
        probe = resolved / ".coblast_write_test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
    except OSError as exc:
        raise ValueError(f"Cannot write to {resolved}: {exc}") from exc
    return resolved


def validate_data_dir(path: str | Path) -> Path:
    """Validate a candidate *data* (write-to) dir and return it resolved.

    Must be writable AND a filesystem BLAST+ can build databases on: a plain file
    write succeeds on exFAT but makeblastdb does not, so this runs the fast
    FAT/exFAT reject and then an authoritative one-shot makeblastdb probe (when
    BLAST+ is available) — a bad choice fails here rather than mid-run.
    """
    resolved = _validate_writable_dir(path)
    require_blast_capable_data_dir(resolved)
    probe_error = makeblastdb_probe_error(resolved)
    if probe_error:
        raise ValueError(
            f"BLAST+ could not build a test database in {resolved}: {probe_error} "
            "Choose a folder on an NTFS drive that BLAST+ can write databases to."
        )
    return resolved


def validate_sra_reads_dir(path: str | Path) -> Path:
    """Validate a candidate SRA reads (pull-from) dir: writable only.

    No makeblastdb check: prefetch/fasterq-dump write plain .sra/.fasta files,
    which work on exFAT/FAT USB drives. The database step targets the data dir.
    """
    return _validate_writable_dir(path)


def _read_settings() -> dict:
    """Return the settings pointer as a dict (empty on missing/corrupt)."""
    try:
        data = json.loads(_settings_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _write_setting(key: str, value: str | None) -> None:
    """Merge one key into the settings pointer without clobbering its siblings.

    The data dir and SRA reads dir share one settings.json, so a whole-file
    rewrite would wipe the other pointer. `value=None` removes the key.
    """
    settings = _settings_path()
    settings.parent.mkdir(parents=True, exist_ok=True)
    data = _read_settings()
    if value is None:
        data.pop(key, None)
    else:
        data[key] = value
    settings.write_text(json.dumps(data), encoding="utf-8")


def _load_pointer(key: str) -> Path | None:
    value = _read_settings().get(key)
    return Path(value).expanduser().resolve() if value else None


def load_saved_data_dir() -> Path | None:
    """Return the user's saved data dir, or None if unset/unreadable/corrupt."""
    return _load_pointer("data_dir")


def save_data_dir(path: str | Path) -> Path:
    """Validate and persist the chosen data dir to the fixed pointer; return it."""
    validated = validate_data_dir(path)
    _write_setting("data_dir", str(validated))
    return validated


def load_saved_sra_reads_dir() -> Path | None:
    """Return the saved SRA reads (pull-from) dir, or None if unset."""
    return _load_pointer("sra_reads_dir")


def save_sra_reads_dir(path: str | Path) -> Path | None:
    """Persist the SRA reads (pull-from) dir; a blank value clears it.

    Cleared means "use the data dir" (reads and databases share one drive). Only
    writability is required, so an exFAT/USB drive is allowed here.
    """
    if not str(path).strip():
        _write_setting("sra_reads_dir", None)
        return None
    validated = validate_sra_reads_dir(path)
    _write_setting("sra_reads_dir", str(validated))
    return validated


def sra_reads_dir() -> Path:
    """Effective SRA reads (pull-from) root where prefetch/fasterq-dump write.

    Precedence: COBLAST_SRA_READS_DIR env, then the saved pull-from pointer (when
    its drive is present), else ``<data_dir>/sra``. Defaulting to the data dir
    means reads and databases share one drive out of the box — the same-drive
    case, which already passed the write-to check. Read live on every fetch and
    scan, so a change here takes effect without a restart.
    """
    env = os.environ.get("COBLAST_SRA_READS_DIR")
    if env:
        return Path(env).expanduser().resolve()
    saved = load_saved_sra_reads_dir()
    if saved and Path(saved.anchor).exists():
        return saved
    return runtime_data_dir() / "sra"


def runtime_data_dir() -> Path:
    """Choose where mutable app data should live."""
    env_data_dir = os.environ.get("COBLAST_DATA_DIR")
    if env_data_dir:
        return Path(env_data_dir).expanduser().resolve()
    saved = load_saved_data_dir()
    if saved and Path(saved.anchor).exists():
        return saved
    # ponytail: a saved dir on a gone drive (unplugged USB, or one that came back
    # as a different letter after reboot) is ignored so the app still launches
    # from the default location instead of dying in setup_data_location's mkdir.
    # Self-heals when the drive returns. Add a picker re-prompt here if support
    # tickets show users need to be asked rather than silently relocated.
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
    bin_dir = ensure_tool_bin("blast", blast_bin_dir)
    exe_path = bin_dir / tool_name(name)
    if not exe_path.exists():
        raise FileNotFoundError(
            f"Could not find {exe_path}. Set BLAST_BIN to your BLAST+ bin directory."
        )
    return exe_path


# CAP3 contig assembler (optional; used by the eToL re-probing workflow). Like
# BLAST+, the binary is resolved from an env override, then a bundled copy, then
# PATH, so the same packaging story works from source and from the frozen .exe.
# Unlike BLAST+, CAP3 is optional: callers gate on availability and skip the
# assembly step when no binary is present rather than failing the whole run.
def _ugene_cap3_candidates() -> list[Path]:
    """Return the ``tools/cap3`` folders of a default UGENE install.

    UGENE ships the CAP3 assembler under ``<install>/tools/cap3`` but does not
    add it to PATH, so probing the standard Program Files locations lets a
    stock UGENE install be found with no ``CAP3_BIN`` setup on the user's part.
    CAP3 is not redistributed with COBLAST+ for licensing reasons; users install
    UGENE (which provides CAP3) instead.
    """
    bases: list[Path] = []
    for env_var in ("ProgramW6432", "ProgramFiles", "ProgramFiles(x86)"):
        value = os.environ.get(env_var)
        if value:
            bases.append(Path(value))
    # Fall back to the conventional install roots if those env vars are unset.
    bases += [Path(r"C:\Program Files"), Path(r"C:\Program Files (x86)")]

    seen: set[Path] = set()
    candidates: list[Path] = []
    for base in bases:
        cap3_dir = base / "Unipro UGENE" / "tools" / "cap3"
        if cap3_dir not in seen:
            seen.add(cap3_dir)
            candidates.append(cap3_dir)
    return candidates


def cap3_bin_dir() -> Path | None:
    """Locate a directory containing the CAP3 binary, or None to fall back to PATH.

    Precedence: the ``CAP3_BIN`` environment variable (a directory, like
    ``BLAST_BIN``), then a bundled ``cap3/bin`` folder inside the resource root,
    then a default Unipro UGENE install (which ships CAP3 under ``tools/cap3``).
    """
    env_cap3_bin = os.environ.get("CAP3_BIN")
    if env_cap3_bin:
        return Path(env_cap3_bin)

    bundled_cap3_bin = resource_path("cap3", "bin")
    if bundled_cap3_bin.exists():
        return bundled_cap3_bin

    # A default UGENE install ships CAP3 but does not put it on PATH; probe its
    # standard tools/cap3 location so users only have to install UGENE.
    for candidate in _ugene_cap3_candidates():
        if (candidate / tool_name("cap3")).exists():
            return candidate

    return None


def cap3_exe() -> Path:
    """Resolve the CAP3 executable from CAP3_BIN, a bundled copy, or PATH.

    Raises ``FileNotFoundError`` when CAP3 cannot be found anywhere, with a hint
    on how to supply it. Callers that treat assembly as optional should check
    availability first rather than letting this propagate.
    """
    directory = cap3_bin_dir()
    if directory is not None:
        candidate = directory / tool_name("cap3")
        if candidate.exists():
            return candidate

    found = shutil.which("cap3") or shutil.which(tool_name("cap3"))
    if found:
        return Path(found)

    raise FileNotFoundError(
        "Could not find the CAP3 assembler. Set CAP3_BIN to the directory "
        "containing the cap3 executable, bundle it under cap3/bin, or add cap3 "
        "to PATH."
    )


# --- Optional auto-install of external tools -------------------------------
#
# When a required binary is not already present (env override, bundled copy, or
# a known install location), COBLAST can fetch a portable build and unpack it
# into the per-user data dir — no admin rights, fully reversible (delete the
# folder). CAP3 is deliberately absent: its license forbids redistribution, so
# it stays detect-only (install UGENE).
#
#   proof:    the executable whose presence proves an unpacked tree is usable.
#   url:      pinned version + host (never "LATEST"); HTTPS is enforced.
#   sha256:   pin for real supply-chain integrity. None here means "verify
#             against NCBI's published .md5 sidecar instead" — that catches a
#             corrupt download, with trust rooted in the same TLS host. A
#             download with neither a pinned hash nor a sidecar fails closed.
#   bin_glob: where `proof` lives under the extracted archive.
#
# ponytail: pin real sha256 values before shipping (defense in depth) and
# confirm the URLs still resolve; NCBI bumps versions.
_DOWNLOADABLE_TOOLS: dict[str, dict[str, str | None]] = {
    "blast": {
        "proof": "blastn",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/"
            "ncbi-blast-2.17.0+-x64-win64.tar.gz"
        ),
        "sha256": None,
        "bin_glob": "ncbi-blast-*/bin",
    },
    "sra": {
        "proof": "fastq-dump",
        "url": (
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.0/"
            "sratoolkit.3.2.0-win64.zip"
        ),
        # NCBI ships no .md5 sidecar for the SRA sdk, so this must be pinned or
        # the download fails closed. Computed from the 3.2.0 win64 zip.
        "sha256": "4b090fc4e12f21203909b997ebdc1140de378d959e4c898f514d6f31d915702e",
        "bin_glob": "sratoolkit.*/bin",
    },
}


def _hash_file(path: Path, algo: str) -> str:
    digest = hashlib.new(algo)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_download(archive: Path, url: str, sha256: str | None) -> None:
    """Fail closed unless the archive matches a pinned sha256 or NCBI's .md5."""
    if sha256:
        actual = _hash_file(archive, "sha256")
        if actual.lower() != sha256.lower():
            archive.unlink(missing_ok=True)
            raise RuntimeError(f"sha256 mismatch for {archive.name} (got {actual}).")
        return
    try:
        with urllib.request.urlopen(url + ".md5", timeout=30) as response:
            expected = response.read().decode("utf-8", "replace").split()[0]
    except Exception as exc:  # no sidecar -> refuse to run an unverified binary
        archive.unlink(missing_ok=True)
        raise RuntimeError(
            f"No pinned sha256 and no .md5 sidecar for {url}; refusing to run an "
            "unverified download. Pin sha256 in _DOWNLOADABLE_TOOLS."
        ) from exc
    actual = _hash_file(archive, "md5")
    if actual.lower() != expected.lower():
        archive.unlink(missing_ok=True)
        raise RuntimeError(f"md5 mismatch for {archive.name} (got {actual}).")


def _extract_archive(archive: Path, dest: Path) -> None:
    name = archive.name.lower()
    if name.endswith(".zip"):
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(dest)
    elif name.endswith((".tar.gz", ".tgz")):
        with tarfile.open(archive, "r:gz") as bundle:
            try:
                bundle.extractall(dest, filter="data")  # 3.12+: block path traversal
            except TypeError:
                bundle.extractall(dest)  # ponytail: <3.12; source is checksum-verified NCBI
    else:
        raise RuntimeError(f"Unknown archive type: {archive.name}")


def _find_tool_bin(root: Path, bin_glob: str, proof_exe: str) -> Path | None:
    if not root.exists():
        return None
    for candidate in sorted(root.glob(bin_glob)):
        if (candidate / proof_exe).exists():
            return candidate
    return None


def download_verify_extract(url: str, dest: Path, sha256: str | None = None) -> None:
    """Download `url` over HTTPS, verify it, then unpack it into `dest`."""
    if not url.lower().startswith("https://"):
        raise ValueError(f"Refusing non-HTTPS download URL: {url}")
    dest.mkdir(parents=True, exist_ok=True)
    archive = dest / url.rsplit("/", 1)[-1]
    # ponytail: stdlib one-shot download; add resume/retry if flaky networks bite.
    urllib.request.urlretrieve(url, archive)
    _verify_download(archive, url, sha256)
    _extract_archive(archive, dest)
    archive.unlink(missing_ok=True)


def ensure_tool_bin(name: str, detector) -> Path | None:
    """Return a usable bin dir for `name`, auto-installing it if possible.

    `detector` is the tool's existing resolver (env -> bundled -> known install
    locations). If it points at a real install, that wins. Otherwise, if the
    tool is downloadable, fetch a portable build into the data dir once and
    reuse it thereafter. Returns None only when the tool is neither installed
    nor downloadable (e.g. CAP3).
    """
    spec = _DOWNLOADABLE_TOOLS.get(name)
    proof_exe = tool_name(spec["proof"]) if spec else None

    found = detector()
    if found is not None and (proof_exe is None or (Path(found) / proof_exe).exists()):
        return Path(found)
    if spec is None:  # not auto-installable (CAP3): hand back whatever detection found
        return Path(found) if found is not None else None

    install_root = runtime_data_dir() / "tools" / name
    bin_dir = _find_tool_bin(install_root, spec["bin_glob"], proof_exe)
    if bin_dir is None:  # not fetched yet
        download_verify_extract(spec["url"], install_root, spec["sha256"])
        bin_dir = _find_tool_bin(install_root, spec["bin_glob"], proof_exe)
    if bin_dir is None:
        raise RuntimeError(f"Installed {name} but no {proof_exe} under {install_root}.")
    return bin_dir


def verify_downloadable_tools() -> None:
    """Pre-ship gate: assert every auto-install entry is fetchable + verifiable.

    Run `python config.py --check-downloads` before handing testers a build (or
    wire it into CI). It catches exactly what a tester would otherwise hit at
    runtime: a URL that 404s, or an entry with neither a pinned sha256 nor a
    .md5 sidecar (which fails closed). HEAD requests only — no large download.
    """
    for name, spec in _DOWNLOADABLE_TOOLS.items():
        url = spec["url"]
        try:
            urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30).close()
        except Exception as exc:
            raise AssertionError(f"{name}: download URL not reachable: {url} ({exc})")
        if spec["sha256"]:
            continue
        try:
            urllib.request.urlopen(url + ".md5", timeout=30).close()
        except Exception:
            raise AssertionError(
                f"{name}: no pinned sha256 and no .md5 sidecar for {url}; testers "
                "would hit a fail-closed download. Pin sha256 in _DOWNLOADABLE_TOOLS."
            )
    print("downloadable tools OK")


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


if __name__ == "__main__":
    # ponytail: offline checks for the auto-install logic (no network touched).
    import tempfile

    # HTTPS is enforced before anything is fetched.
    try:
        download_verify_extract("http://example.com/x.zip", Path(tempfile.gettempdir()))
    except ValueError:
        pass
    else:
        raise AssertionError("non-HTTPS URL should be refused")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        assert _find_tool_bin(root, "sratoolkit.*/bin", "fastq-dump") is None
        bin_dir = root / "sratoolkit.3.2.0" / "bin"
        bin_dir.mkdir(parents=True)
        (bin_dir / "fastq-dump").write_text("x")
        assert _find_tool_bin(root, "sratoolkit.*/bin", "fastq-dump") == bin_dir

        # A pinned sha256 that does not match must fail closed.
        blob = root / "blob.zip"
        blob.write_bytes(b"hello")
        try:
            _verify_download(blob, "https://x/blob.zip", "deadbeef")
        except RuntimeError:
            pass
        else:
            raise AssertionError("sha256 mismatch should fail closed")

    # CAP3 is not auto-installable: ensure_tool_bin returns detection, never downloads.
    assert ensure_tool_bin("cap3", lambda: None) is None
    assert ensure_tool_bin("cap3", lambda: Path("/opt/ugene")) == Path("/opt/ugene")
    print("config auto-install checks OK")

    # ponytail: data-dir pointer + validation checks (no network, no real data dir).
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        saved_env = {k: os.environ.get(k) for k in ("LOCALAPPDATA", "XDG_DATA_HOME")}
        os.environ["LOCALAPPDATA"] = str(base)  # Windows pointer home
        os.environ["XDG_DATA_HOME"] = str(base)  # Linux pointer home
        try:
            assert load_saved_data_dir() is None
            # load reads back a hand-written pointer (read path does not validate)
            _settings_path().parent.mkdir(parents=True, exist_ok=True)
            _settings_path().write_text(
                json.dumps({"data_dir": str(base / "store")}), encoding="utf-8"
            )
            assert load_saved_data_dir() == (base / "store").expanduser().resolve()
            # a corrupt pointer degrades to None instead of crashing
            _settings_path().write_text("{ not json", encoding="utf-8")
            assert load_saved_data_dir() is None
            # spaces are always rejected, before any disk is touched
            try:
                validate_data_dir("some folder/with space")
            except ValueError:
                pass
            else:
                raise AssertionError("spaced data dir should be rejected")
            # a clean, writable dir round-trips through save/load
            if " " not in str(base):
                chosen = save_data_dir(base / "picked")
                assert chosen == (base / "picked").resolve()
                assert load_saved_data_dir() == chosen
                # the SRA reads pointer coexists with data_dir (merge, no clobber)
                reads = save_sra_reads_dir(base / "reads")
                assert reads == (base / "reads").resolve()
                assert load_saved_sra_reads_dir() == reads
                assert load_saved_data_dir() == chosen  # data_dir survived the reads write
                # clearing the reads pointer leaves data_dir intact
                assert save_sra_reads_dir("") is None
                assert load_saved_sra_reads_dir() is None
                assert load_saved_data_dir() == chosen
            # the FAT/exFAT check must not false-positive on a normal (NTFS) temp dir
            assert blast_incompatible_filesystem(base) is None
            try:
                require_blast_capable_data_dir(base)
            except ValueError:
                raise AssertionError("NTFS temp dir wrongly rejected as BLAST-incompatible")
            # a saved pointer on a missing drive (unplugged USB) is ignored, not returned
            if os.name == "nt":
                for letter in "QZYXWVUT":
                    if not Path(f"{letter}:\\").exists():
                        _settings_path().write_text(
                            json.dumps({"data_dir": f"{letter}:\\coblast"}), encoding="utf-8"
                        )
                        assert runtime_data_dir() != Path(f"{letter}:\\coblast")
                        break
        finally:
            for key, value in saved_env.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
    print("data-dir pointer checks OK")

    # Opt-in network check (run before shipping a test build): python config.py --check-downloads
    if "--check-downloads" in sys.argv:
        verify_downloadable_tools()
