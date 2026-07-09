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


def insert_database(
    *,
    display_name: str,
    prefix: Path,
    category: str = "custom",
    status: str = "available",
):
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
            VALUES (?, 'nucl', ?, ?, ?, ?, ?)
            """,
            (display_name, normalized_prefix, category, now, now, status),
        )
        database_id = cursor.lastrowid
    return registry.get_database(database_id)


def test_create_database_passes_spaceless_in_arg_for_spaced_source_dir(
    isolated_registry, monkeypatch
):
    # A FASTA inside a directory with a space must never reach makeblastdb's -in as
    # a spaced string: BLAST+ splits -in on whitespace, so "E:\S DRIVE\reads.fasta"
    # would be read as two files. We run from the dir and pass the bare filename.
    source_dir = isolated_registry / "S DRIVE"
    source_dir.mkdir()
    source = source_dir / "reads.fasta"
    source.write_text(">r\nACGT\n", encoding="ascii")

    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["cwd"] = kwargs.get("cwd")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(registry, "blast_exe", lambda name: Path(name))
    monkeypatch.setattr(registry.subprocess, "run", fake_run)
    monkeypatch.setattr(registry, "upsert_database", lambda **kwargs: SimpleNamespace(**kwargs))

    registry.create_database_from_fasta(
        display_name="Spaced Source",
        db_type="nucl",
        source_fasta_path=source,
    )

    in_value = captured["command"][captured["command"].index("-in") + 1]
    assert in_value == "reads.fasta"
    assert " " not in in_value
    assert Path(captured["cwd"]) == source_dir


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


def test_remove_missing_databases_removes_only_missing_rows(isolated_registry):
    missing_one = insert_database(
        display_name="Missing A",
        prefix=isolated_registry / "missing-a",
        status="missing",
    )
    available = insert_database(
        display_name="Available",
        prefix=isolated_registry / "available",
    )
    missing_two = insert_database(
        display_name="Missing B",
        prefix=isolated_registry / "missing-b",
        status="missing",
    )
    invalid = insert_database(
        display_name="Invalid",
        prefix=isolated_registry / "invalid",
        status="invalid",
    )

    removed = registry.remove_missing_databases()

    assert [database.display_name for database in removed] == ["Missing A", "Missing B"]
    assert [database.id for database in registry.list_databases()] == [
        available.id,
        invalid.id,
    ]
    assert registry.database_prefix_was_removed(missing_one.db_prefix_path)
    assert registry.database_prefix_was_removed(missing_two.db_prefix_path)
    assert not registry.database_prefix_was_removed(available.db_prefix_path)
    assert not registry.database_prefix_was_removed(invalid.db_prefix_path)


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


def test_remove_missing_route_reports_removed_count(monkeypatch):
    monkeypatch.setattr(
        app_module,
        "remove_missing_databases",
        lambda: [SimpleNamespace(), SimpleNamespace()],
    )
    client = app_module.app.test_client()

    response = client.post("/databases/remove-missing")

    assert response.status_code == 302
    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["message"] == [
        "Removed 2 missing databases from the registry. BLAST files were not deleted."
    ]


def test_remove_missing_route_handles_empty_selection(monkeypatch):
    monkeypatch.setattr(app_module, "remove_missing_databases", lambda: [])
    client = app_module.app.test_client()

    response = client.post("/databases/remove-missing")

    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["message"] == ["No missing databases were found."]


def test_verify_selected_route_verifies_each_id(monkeypatch):
    verified = []
    monkeypatch.setattr(
        app_module,
        "verify_database",
        lambda database_id: verified.append(database_id)
        or SimpleNamespace(display_name=f"DB {database_id}", status="available"),
    )
    client = app_module.app.test_client()

    response = client.post(
        "/databases/verify-selected", data={"selected_db": ["3", "5"]}
    )

    assert response.status_code == 302
    assert verified == [3, 5]
    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["message"] == ["Verified 2 selected databases."]


def test_verify_selected_route_requires_a_selection():
    client = app_module.app.test_client()

    response = client.post("/databases/verify-selected", data={})

    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["error"] == ["Select at least one database to verify."]


def test_remove_selected_route_removes_each_id(monkeypatch):
    removed = []
    monkeypatch.setattr(
        app_module,
        "remove_database",
        lambda database_id: removed.append(database_id)
        or SimpleNamespace(display_name=f"DB {database_id}"),
    )
    client = app_module.app.test_client()

    response = client.post(
        "/databases/remove-selected", data={"selected_db": ["7"]}
    )

    assert response.status_code == 302
    assert removed == [7]
    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["message"] == [
        "Removed 1 selected database from the registry. BLAST files were not deleted."
    ]


def test_remove_selected_route_requires_a_selection():
    client = app_module.app.test_client()

    response = client.post("/databases/remove-selected", data={})

    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["error"] == ["Select at least one database to remove."]


def test_remove_selected_route_reports_per_database_errors(monkeypatch):
    def fake_remove(database_id):
        if database_id == 2:
            raise ValueError("not found")
        return SimpleNamespace(display_name=f"DB {database_id}")

    monkeypatch.setattr(app_module, "remove_database", fake_remove)
    client = app_module.app.test_client()

    response = client.post(
        "/databases/remove-selected", data={"selected_db": ["1", "2"]}
    )

    query = parse_qs(urlparse(response.headers["Location"]).query)
    assert query["message"] == [
        "Removed 1 selected database from the registry. BLAST files were not deleted."
    ]
    assert query["error"] == ["Database 2: not found"]


def test_database_page_includes_remove_confirmation(monkeypatch):
    database = SimpleNamespace(
        id=7,
        display_name="Clinical DB",
        description="",
        db_type="nucl",
        category="custom",
        status="missing",
        sequence_count=10,
        last_verified_at="2026-06-13T12:00:00+00:00",
        db_prefix_path=r"C:\BLAST\clinical",
    )
    monkeypatch.setattr(app_module, "ensure_demo_databases", lambda: None)
    monkeypatch.setattr(app_module, "list_databases", lambda: [database])
    client = app_module.app.test_client()

    response = client.get("/databases")

    assert response.status_code == 200
    # Each row exposes a checkbox carrying the database id for bulk selection.
    assert b'name="selected_db"' in response.data
    assert b'value="7"' in response.data
    # Bulk verify/remove buttons replace the old per-row actions.
    assert b"Verify Selected" in response.data
    assert b"Remove Selected" in response.data
    assert b'data-remove-missing-count="1"' in response.data
    assert b"Remove All Missing (1)" in response.data
    assert b"The FASTA and BLAST database files will remain on disk." in response.data
