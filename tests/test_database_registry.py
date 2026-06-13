"""Tests for removing databases from the local SQLite registry."""

from pathlib import Path
from types import SimpleNamespace
from urllib.parse import parse_qs, urlparse

import pytest

import app as app_module
import database_registry as registry


@pytest.fixture
def isolated_registry(tmp_path, monkeypatch):
    instance_dir = tmp_path / "instance"
    sample_data_dir = tmp_path / "sample_data"
    sample_data_dir.mkdir()
    (sample_data_dir / "toy_nt.fasta").write_text(">toy_nt\nACGT\n", encoding="ascii")
    (sample_data_dir / "toy_protein.fasta").write_text(">toy_protein\nMKT\n", encoding="ascii")

    monkeypatch.setattr(registry, "INSTANCE_DIR", instance_dir)
    monkeypatch.setattr(registry, "REGISTRY_PATH", instance_dir / "registry.sqlite")
    monkeypatch.setattr(registry, "MANAGED_DATABASE_DIR", instance_dir / "databases")
    monkeypatch.setattr(registry, "SAMPLE_DATA_DIR", sample_data_dir)
    registry.init_registry()
    return instance_dir


def insert_database(*, display_name: str, prefix: Path, category: str = "custom"):
    normalized_prefix = registry.blast_safe_path(prefix)
    now = registry.utc_now()
    with registry.registry_connection() as conn:
        cursor = conn.execute(
            """
            INSERT INTO blast_databases (
                display_name,
                db_type,
                db_prefix_path,
                category,
                created_at,
                last_verified_at,
                status
            )
            VALUES (?, 'nucl', ?, ?, ?, ?, 'available')
            """,
            (display_name, normalized_prefix, category, now, now),
        )
        database_id = cursor.lastrowid
    return registry.get_database(database_id)


def test_remove_database_deletes_only_the_registry_row(isolated_registry):
    prefix = isolated_registry / "external" / "patient"
    prefix.parent.mkdir()
    blast_index = prefix.with_suffix(".nin")
    blast_index.write_text("index data", encoding="ascii")
    database = insert_database(display_name="Patient DB", prefix=prefix)

    removed = registry.remove_database(database.id)

    assert removed.display_name == "Patient DB"
    assert registry.list_databases() == []
    assert registry.database_prefix_was_removed(prefix)
    assert blast_index.read_text(encoding="ascii") == "index data"


def test_remove_database_rejects_unknown_id(isolated_registry):
    with pytest.raises(ValueError, match="No registered database exists with ID 999"):
        registry.remove_database(999)


def test_removed_demo_database_is_not_seeded_again(isolated_registry, monkeypatch):
    toy_prefix = registry.MANAGED_DATABASE_DIR / "toy_nt"
    database = insert_database(
        display_name="Toy Nucleotide Test Database",
        prefix=toy_prefix,
        category="toy",
    )
    registry.remove_database(database.id)
    created_prefixes = []

    def record_creation(**database_fields):
        created_prefixes.append(Path(database_fields["db_prefix_path"]).name)

    monkeypatch.setattr(registry, "create_database_from_fasta", record_creation)

    registry.ensure_demo_databases()

    assert "toy_nt" not in created_prefixes
    assert "toy_protein" in created_prefixes


def test_explicit_registration_restores_a_removed_prefix(isolated_registry, monkeypatch):
    prefix = isolated_registry / "external" / "viral"
    database = insert_database(display_name="Old viral DB", prefix=prefix)
    registry.remove_database(database.id)
    monkeypatch.setattr(
        registry,
        "verify_database_prefix",
        lambda _prefix: {
            "status": "available",
            "database_title": "Viral panel",
            "sequence_count": 12,
            "notes": "",
        },
    )
    monkeypatch.setattr(registry, "blast_version", lambda: "test-version")

    restored = registry.register_existing_database(
        display_name="Restored viral DB",
        db_type="nucl",
        db_prefix_path=prefix,
    )

    assert restored.display_name == "Restored viral DB"
    assert not registry.database_prefix_was_removed(prefix)


def test_remove_route_names_database_and_explains_files_remain(monkeypatch):
    monkeypatch.setattr(
        app_module,
        "remove_database",
        lambda database_id: SimpleNamespace(display_name=f"Patient DB {database_id}"),
    )
    client = app_module.app.test_client()

    response = client.post("/databases/42/remove")

    assert response.status_code == 302
    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["message"] == [
        "Patient DB 42 was removed from the registry. BLAST files were not deleted."
    ]


def test_database_page_includes_remove_confirmation(monkeypatch):
    database = SimpleNamespace(
        id=7,
        display_name="Clinical DB",
        description="",
        db_type="nucl",
        category="custom",
        status="available",
        sequence_count=10,
        last_verified_at="2026-06-13T12:00:00+00:00",
        db_prefix_path=r"C:\BLAST\clinical",
    )
    monkeypatch.setattr(app_module, "ensure_demo_databases", lambda: None)
    monkeypatch.setattr(app_module, "list_databases", lambda: [database])
    client = app_module.app.test_client()

    response = client.get("/databases")

    assert response.status_code == 200
    assert b'data-remove-database="Clinical DB"' in response.data
    assert b"The FASTA and BLAST database files will remain on disk." in response.data
