from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import ctypes
import os
from pathlib import Path
import re
import shlex
import sqlite3
import subprocess

from config import blast_exe


PROJECT_ROOT = Path(__file__).resolve().parent
INSTANCE_DIR = PROJECT_ROOT / "instance"
REGISTRY_PATH = INSTANCE_DIR / "database_registry.sqlite"
MANAGED_DATABASE_DIR = INSTANCE_DIR / "databases"
SAMPLE_DATA_DIR = PROJECT_ROOT / "sample_data"

DB_TYPES = {"nucl", "prot"}
DB_CATEGORIES = {"viral", "human", "eToL-V", "toy", "custom"}
DB_STATUSES = {"available", "missing", "invalid"}


@dataclass(frozen=True)
class RegisteredDatabase:
    id: int
    display_name: str
    db_type: str
    db_prefix_path: str
    source_fasta_path: str
    description: str
    category: str
    created_at: str
    last_verified_at: str
    blast_version: str
    makeblastdb_command: str
    sequence_count: int | None
    database_title: str
    status: str
    notes: str


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def registry_connection() -> sqlite3.Connection:
    INSTANCE_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(REGISTRY_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_registry() -> None:
    with registry_connection() as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS blast_databases (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                display_name TEXT NOT NULL,
                db_type TEXT NOT NULL CHECK (db_type IN ('nucl', 'prot')),
                db_prefix_path TEXT NOT NULL UNIQUE,
                source_fasta_path TEXT NOT NULL DEFAULT '',
                description TEXT NOT NULL DEFAULT '',
                category TEXT NOT NULL DEFAULT 'custom',
                created_at TEXT NOT NULL,
                last_verified_at TEXT NOT NULL DEFAULT '',
                blast_version TEXT NOT NULL DEFAULT '',
                makeblastdb_command TEXT NOT NULL DEFAULT '',
                sequence_count INTEGER,
                database_title TEXT NOT NULL DEFAULT '',
                status TEXT NOT NULL CHECK (status IN ('available', 'missing', 'invalid')),
                notes TEXT NOT NULL DEFAULT ''
            )
            """
        )


def row_to_database(row: sqlite3.Row) -> RegisteredDatabase:
    return RegisteredDatabase(
        id=row["id"],
        display_name=row["display_name"],
        db_type=row["db_type"],
        db_prefix_path=row["db_prefix_path"],
        source_fasta_path=row["source_fasta_path"],
        description=row["description"],
        category=row["category"],
        created_at=row["created_at"],
        last_verified_at=row["last_verified_at"],
        blast_version=row["blast_version"],
        makeblastdb_command=row["makeblastdb_command"],
        sequence_count=row["sequence_count"],
        database_title=row["database_title"],
        status=row["status"],
        notes=row["notes"],
    )


def validate_db_type(db_type: str) -> str:
    cleaned = db_type.strip().lower()
    if cleaned not in DB_TYPES:
        raise ValueError("Database type must be 'nucl' or 'prot'.")
    return cleaned


def validate_category(category: str) -> str:
    cleaned = category.strip() or "custom"
    if cleaned not in DB_CATEGORIES:
        raise ValueError(
            "Database category must be one of: " + ", ".join(sorted(DB_CATEGORIES))
        )
    return cleaned


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", value.strip()).strip("_").lower()
    return slug or "blast_database"


def get_windows_short_path(path: Path) -> str:
    if os.name != "nt":
        return str(path)

    try:
        get_short_path_name = ctypes.windll.kernel32.GetShortPathNameW
    except AttributeError:
        return str(path)

    source = str(path)
    length = get_short_path_name(source, None, 0)
    if length == 0:
        return source

    buffer = ctypes.create_unicode_buffer(length)
    if get_short_path_name(source, buffer, length) == 0:
        return source
    return buffer.value


def blast_safe_path(path: str | Path) -> str:
    resolved = Path(path).expanduser().resolve()
    if os.name != "nt" or " " not in str(resolved):
        return str(resolved)

    if resolved.exists():
        return get_windows_short_path(resolved)

    parent = resolved.parent
    if parent.exists():
        return str(Path(get_windows_short_path(parent)) / resolved.name)

    return str(resolved)


def command_to_string(command: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def blast_version() -> str:
    completed = subprocess.run(
        [str(blast_exe("blastn")), "-version"],
        capture_output=True,
        text=True,
        check=False,
    )
    return (completed.stdout or completed.stderr).strip().replace("\n", " ")


def parse_blastdbcmd_info(stdout: str) -> tuple[str, int | None]:
    title = ""
    sequence_count: int | None = None

    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if line.startswith("Database:"):
            title = line.removeprefix("Database:").strip()
            continue
        match = re.search(r"([0-9,]+)\s+sequences?;", line)
        if match:
            sequence_count = int(match.group(1).replace(",", ""))

    return title, sequence_count


def verify_database_prefix(db_prefix_path: str | Path) -> dict[str, str | int | None]:
    db_prefix = blast_safe_path(db_prefix_path)
    completed = subprocess.run(
        [str(blast_exe("blastdbcmd")), "-db", db_prefix, "-info"],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        message = (completed.stderr or completed.stdout).strip()
        status = "missing" if "No alias or index file found" in message else "invalid"
        return {
            "status": status,
            "database_title": "",
            "sequence_count": None,
            "notes": message,
        }

    title, sequence_count = parse_blastdbcmd_info(completed.stdout)
    return {
        "status": "available",
        "database_title": title,
        "sequence_count": sequence_count,
        "notes": "",
    }


def list_databases() -> list[RegisteredDatabase]:
    init_registry()
    with registry_connection() as conn:
        rows = conn.execute(
            "SELECT * FROM blast_databases ORDER BY display_name COLLATE NOCASE"
        ).fetchall()
    return [row_to_database(row) for row in rows]


def list_compatible_databases(required_db_type: str) -> list[RegisteredDatabase]:
    required_db_type = validate_db_type(required_db_type)
    return [
        database
        for database in list_databases()
        if database.db_type == required_db_type and database.status == "available"
    ]


def get_database_by_prefix(db_prefix_path: str | Path) -> RegisteredDatabase | None:
    init_registry()
    prefix = blast_safe_path(db_prefix_path)
    with registry_connection() as conn:
        row = conn.execute(
            "SELECT * FROM blast_databases WHERE db_prefix_path = ?",
            (prefix,),
        ).fetchone()
    return row_to_database(row) if row else None


def get_database(database_id: int) -> RegisteredDatabase:
    init_registry()
    with registry_connection() as conn:
        row = conn.execute(
            "SELECT * FROM blast_databases WHERE id = ?",
            (database_id,),
        ).fetchone()
    if row is None:
        raise ValueError(f"No registered database exists with ID {database_id}.")
    return row_to_database(row)


def upsert_database(
    *,
    display_name: str,
    db_type: str,
    db_prefix_path: str | Path,
    source_fasta_path: str | Path | None = None,
    description: str = "",
    category: str = "custom",
    makeblastdb_command: str = "",
    notes: str = "",
) -> RegisteredDatabase:
    init_registry()
    if not display_name.strip():
        raise ValueError("Enter a database display name.")
    if not str(db_prefix_path).strip():
        raise ValueError("Enter a BLAST database prefix path.")

    db_type = validate_db_type(db_type)
    category = validate_category(category)
    prefix = blast_safe_path(db_prefix_path)
    source = str(Path(source_fasta_path).expanduser().resolve()) if source_fasta_path else ""
    verified = verify_database_prefix(prefix)
    now = utc_now()

    status = str(verified["status"])
    registry_notes = notes.strip()
    if verified["notes"]:
        registry_notes = "\n".join(part for part in [registry_notes, str(verified["notes"])] if part)

    with registry_connection() as conn:
        conn.execute(
            """
            INSERT INTO blast_databases (
                display_name,
                db_type,
                db_prefix_path,
                source_fasta_path,
                description,
                category,
                created_at,
                last_verified_at,
                blast_version,
                makeblastdb_command,
                sequence_count,
                database_title,
                status,
                notes
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(db_prefix_path) DO UPDATE SET
                display_name = excluded.display_name,
                db_type = excluded.db_type,
                source_fasta_path = excluded.source_fasta_path,
                description = excluded.description,
                category = excluded.category,
                last_verified_at = excluded.last_verified_at,
                blast_version = excluded.blast_version,
                makeblastdb_command = excluded.makeblastdb_command,
                sequence_count = excluded.sequence_count,
                database_title = excluded.database_title,
                status = excluded.status,
                notes = excluded.notes
            """,
            (
                display_name.strip(),
                db_type,
                prefix,
                source,
                description.strip(),
                category,
                now,
                now,
                blast_version(),
                makeblastdb_command,
                verified["sequence_count"],
                verified["database_title"],
                status,
                registry_notes,
            ),
        )
        database_id = conn.execute(
            "SELECT id FROM blast_databases WHERE db_prefix_path = ?",
            (prefix,),
        ).fetchone()["id"]
    return get_database(database_id)


def register_existing_database(
    *,
    display_name: str,
    db_type: str,
    db_prefix_path: str | Path,
    source_fasta_path: str | Path | None = None,
    description: str = "",
    category: str = "custom",
    notes: str = "",
) -> RegisteredDatabase:
    return upsert_database(
        display_name=display_name,
        db_type=db_type,
        db_prefix_path=db_prefix_path,
        source_fasta_path=source_fasta_path,
        description=description,
        category=category,
        notes=notes,
    )


def create_database_from_fasta(
    *,
    display_name: str,
    db_type: str,
    source_fasta_path: str | Path,
    db_prefix_path: str | Path | None = None,
    description: str = "",
    category: str = "custom",
    notes: str = "",
) -> RegisteredDatabase:
    db_type = validate_db_type(db_type)
    category = validate_category(category)
    if not display_name.strip():
        raise ValueError("Enter a database display name.")
    if not str(source_fasta_path).strip():
        raise ValueError("Enter the source FASTA path.")

    source = Path(source_fasta_path).expanduser().resolve()
    if not source.exists():
        raise ValueError(f"Source FASTA does not exist: {source}")

    if db_prefix_path:
        prefix = Path(db_prefix_path).expanduser().resolve()
    else:
        prefix = MANAGED_DATABASE_DIR / slugify(display_name)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    command = [
        str(blast_exe("makeblastdb")),
        "-in",
        blast_safe_path(source),
        "-dbtype",
        db_type,
        "-out",
        blast_safe_path(prefix),
    ]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError((completed.stderr or completed.stdout).strip())

    return upsert_database(
        display_name=display_name,
        db_type=db_type,
        db_prefix_path=blast_safe_path(prefix),
        source_fasta_path=source,
        description=description,
        category=category,
        makeblastdb_command=command_to_string(command),
        notes=notes,
    )


def verify_database(database_id: int) -> RegisteredDatabase:
    database = get_database(database_id)
    verified = verify_database_prefix(database.db_prefix_path)
    now = utc_now()
    notes = database.notes
    if verified["notes"]:
        notes = str(verified["notes"])

    with registry_connection() as conn:
        conn.execute(
            """
            UPDATE blast_databases
            SET last_verified_at = ?,
                sequence_count = ?,
                database_title = ?,
                status = ?,
                notes = ?
            WHERE id = ?
            """,
            (
                now,
                verified["sequence_count"],
                verified["database_title"],
                verified["status"],
                notes,
                database_id,
            ),
        )
    return get_database(database_id)


def remove_database(database_id: int) -> None:
    init_registry()
    with registry_connection() as conn:
        conn.execute("DELETE FROM blast_databases WHERE id = ?", (database_id,))


def ensure_demo_databases() -> None:
    demo_databases = [
        {
            "display_name": "Toy Nucleotide Test Database",
            "db_type": "nucl",
            "source_fasta_path": SAMPLE_DATA_DIR / "toy_nt.fasta",
            "db_prefix_path": MANAGED_DATABASE_DIR / "toy_nt",
            "description": "Small nucleotide database for testing blastn and tblastn.",
            "category": "toy",
            "notes": "Seeded automatically for local prototype testing.",
        },
        {
            "display_name": "Toy Protein Test Database",
            "db_type": "prot",
            "source_fasta_path": SAMPLE_DATA_DIR / "toy_protein.fasta",
            "db_prefix_path": MANAGED_DATABASE_DIR / "toy_protein",
            "description": "Small protein database for testing blastp and blastx.",
            "category": "toy",
            "notes": "Seeded automatically for local prototype testing.",
        },
    ]

    for database in demo_databases:
        source = Path(database["source_fasta_path"])
        prefix = Path(database["db_prefix_path"])
        if not source.exists():
            continue

        existing = get_database_by_prefix(prefix)
        if existing is not None:
            verified = verify_database(existing.id)
            if verified.status == "available":
                continue

        create_database_from_fasta(**database)
