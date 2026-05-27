"""Flask entry point for the local COBLAST+ web interface.

This module keeps the web layer deliberately thin: routes collect form input,
delegate BLAST/database work to helper modules, and render templates with the
objects those helpers return.
"""

from dataclasses import replace
from pathlib import Path

from flask import Flask, Response, abort, redirect, render_template, request, url_for

from blast_runner import BLAST_PROGRAMS, SENSITIVITY_PRESETS, run_blast
from database_registry import (
    DB_CATEGORIES,
    DB_TYPES,
    create_database_from_fasta,
    ensure_demo_databases,
    get_database,
    list_databases,
    register_existing_database,
    remove_database,
    verify_database,
)
from result_store import (
    batch_rows_as_delimited,
    load_batch_result,
    load_result,
    result_rows_as_delimited,
    save_batch_result,
    save_result,
)
from config import FLASK_HOST, flask_port, resource_path, resource_root
from runtime_estimator import database_storage_bytes, format_bytes, format_duration
from sra_workflow import (
    configured_sra_roots,
    convert_sra_to_pilot_fasta,
    create_pilot_database_from_fasta,
    discover_sra_projects,
    register_sra_blast_database,
    sra_toolkit_bin,
)


app = Flask(
    __name__,
    # resource_path works both from source and from a PyInstaller bundle.
    template_folder=str(resource_path("templates")),
    static_folder=str(resource_path("static")),
)

APOE_PROBE_FASTA = """>AE4_E4=C
CGGACATGGAGGACGTGCGCGGCCGCCTGGTGCAGT
>AE4_E23=T
CGGACATGGAGGACGTGTGCGGCCGCCTGGTGCAGT
>AE2_E34=C
CCGATGACCTGCAGAAGCGCCTGGCAGTGTACCAGG
>AE2_E2=T
CCGATGACCTGCAGAAGTGCCTGGCAGTGTACCAGG
"""
APOE_EXACT_MATCH_FILTER = "100% identity and 100% query coverage"


def numeric_hit_value(hit: dict[str, str], key: str) -> float | None:
    """Parse a numeric hit field that may be blank or formatted text."""
    try:
        return float(hit.get(key, ""))
    except (TypeError, ValueError):
        return None


def filter_apoe_exact_probe_hits(hits: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep only exact APOE probe matches after BLAST output parsing."""
    exact_hits = []
    for hit in hits:
        percent_identity = numeric_hit_value(hit, "pident")
        query_coverage = numeric_hit_value(hit, "qcovs")
        if percent_identity == 100.0 and query_coverage == 100.0:
            exact_hits.append(hit)
    return exact_hits

PROJECT_ROOT = resource_root()


def redirect_to_databases(message: str = "", error: str = ""):
    """Send users back to the database page with optional status text."""
    params = {}
    if message:
        params["message"] = message
    if error:
        params["error"] = error
    return redirect(url_for("databases_page", **params))


def redirect_to_sra(message: str = "", error: str = ""):
    """Send users back to the SRA workbench with optional status text."""
    params = {}
    if message:
        params["message"] = message
    if error:
        params["error"] = error
    return redirect(url_for("sra_page", **params))


def database_options():
    """Attach local storage-size metadata to registered database rows."""
    options = []
    for database in list_databases():
        storage_bytes = database_storage_bytes(database.db_prefix_path)
        options.append(
            {
                "database": database,
                "storage_bytes": storage_bytes,
                "storage_label": format_bytes(storage_bytes),
            }
        )
    return options


@app.get("/")
def index():
    """Render the main search form."""
    try:
        # Demo databases make first-run testing possible without manual setup.
        ensure_demo_databases()
        databases = list_databases()
        registry_error = None
    except Exception as exc:
        databases = []
        registry_error = str(exc)

    return render_template(
        "index.html",
        blast_programs=BLAST_PROGRAMS,
        sensitivity_presets=SENSITIVITY_PRESETS,
        databases=databases,
        registry_error=registry_error,
    )


@app.post("/run-blast")
def run_blast_route():
    """Validate search options, run BLAST locally, and show parsed results."""
    sequence = request.form.get("sequence", "")
    uploaded_query = request.files.get("sequence_file")
    if uploaded_query and uploaded_query.filename:
        # UTF-8 with BOM support keeps FASTA uploads from Windows editors usable.
        sequence = uploaded_query.read().decode("utf-8-sig")

    # Form fields arrive as strings; blast_runner owns numeric validation.
    program = request.form.get("program", "blastn")
    database_id = request.form.get("database_id", "")
    output_format = request.form.get("output_format", "tabular")
    sensitivity_preset = request.form.get("sensitivity_preset", "standard")

    try:
        if program not in BLAST_PROGRAMS:
            raise ValueError(f"Unsupported BLAST program: {program}")
        if not database_id:
            raise ValueError("Choose a registered BLAST database.")

        database = get_database(int(database_id))
        required_db_type = str(BLAST_PROGRAMS[program]["db_type"])
        # The UI filters incompatible databases, but the server still enforces it.
        if database.db_type != required_db_type:
            raise ValueError(
                f"{BLAST_PROGRAMS[program]['label']} requires a {required_db_type} "
                f"database, but {database.display_name} is registered as {database.db_type}."
            )
        if database.status != "available":
            raise ValueError(
                f"{database.display_name} is currently marked as {database.status}. "
                "Verify it on the Databases page before running BLAST."
            )

        # run_blast writes a temporary query file, calls the local BLAST+ binary,
        # parses stdout, and returns a serializable BlastResult dataclass.
        result = run_blast(
            sequence=sequence,
            database=database.db_prefix_path,
            program=program,
            output_format=output_format,
            sensitivity_preset=sensitivity_preset,
            task=request.form.get("task") or None,
            evalue=request.form.get("evalue") or None,
            max_target_seqs=request.form.get("max_target_seqs") or None,
            word_size=request.form.get("word_size") or None,
            perc_identity=request.form.get("perc_identity") or None,
            timeout_seconds=request.form.get("timeout_seconds") or None,
        )
    except Exception as exc:
        return render_template(
            "results.html",
            error=str(exc),
            result=None,
            format_bytes=format_bytes,
            format_duration=format_duration,
        ), 400

    # Persist a JSON copy so the results page can offer CSV/TSV downloads.
    run_id = save_result(result)
    return render_template(
        "results.html",
        error=None,
        result=result,
        run_id=run_id,
        format_bytes=format_bytes,
        format_duration=format_duration,
    )


@app.get("/results/<run_id>.<file_format>")
def download_results(run_id: str, file_format: str):
    """Download a saved result table as CSV or TSV."""
    if file_format not in {"csv", "tsv"}:
        abort(404)

    try:
        result_data = load_result(run_id)
    except FileNotFoundError:
        abort(404)

    delimiter = "," if file_format == "csv" else "\t"
    body = result_rows_as_delimited(result_data, delimiter=delimiter)
    mimetype = "text/csv" if file_format == "csv" else "text/tab-separated-values"
    filename = f"blast_results_{run_id}.{file_format}"
    return Response(
        body,
        mimetype=mimetype,
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )


@app.get("/batch-results/<batch_id>.<file_format>")
def download_batch_results(batch_id: str, file_format: str):
    """Download saved batch results as CSV or TSV."""
    if file_format not in {"csv", "tsv"}:
        abort(404)

    try:
        batch_data = load_batch_result(batch_id)
    except FileNotFoundError:
        abort(404)

    delimiter = "," if file_format == "csv" else "\t"
    body = batch_rows_as_delimited(batch_data, delimiter=delimiter)
    mimetype = "text/csv" if file_format == "csv" else "text/tab-separated-values"
    filename = f"batch_blast_results_{batch_id}.{file_format}"
    return Response(
        body,
        mimetype=mimetype,
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )


@app.get("/batch-blast")
def batch_blast_page():
    """Render the multi-database BLAST form."""
    try:
        ensure_demo_databases()
        options = database_options()
        registry_error = None
    except Exception as exc:
        options = []
        registry_error = str(exc)

    return render_template(
        "batch.html",
        blast_programs=BLAST_PROGRAMS,
        sensitivity_presets=SENSITIVITY_PRESETS,
        database_options=options,
        apoe_probe_fasta=APOE_PROBE_FASTA,
        error=request.args.get("error") or registry_error,
        message=request.args.get("message", ""),
    )


@app.post("/batch-blast")
def run_batch_blast_route():
    """Run one query against multiple registered local databases sequentially."""
    sequence = request.form.get("sequence", "")
    uploaded_query = request.files.get("sequence_file")
    if uploaded_query and uploaded_query.filename:
        sequence = uploaded_query.read().decode("utf-8-sig")

    program = request.form.get("program", "blastn")
    database_ids = request.form.getlist("database_ids")
    output_format = request.form.get("output_format", "tabular")
    sensitivity_preset = request.form.get("sensitivity_preset", "standard")
    apoe_probe_preset = request.form.get("apoe_probe_preset") == "1"
    if apoe_probe_preset:
        sequence = APOE_PROBE_FASTA
        program = "blastn"
        output_format = "tabular"

    if not database_ids:
        return render_template(
            "batch.html",
            blast_programs=BLAST_PROGRAMS,
            sensitivity_presets=SENSITIVITY_PRESETS,
            database_options=database_options(),
            apoe_probe_fasta=APOE_PROBE_FASTA,
            error="Choose at least one database for the batch run.",
            message="",
        ), 400

    database_results = []
    total_runtime_seconds = 0.0
    total_hits = 0
    query_count = 0
    query_total_length = 0

    for raw_database_id in database_ids:
        try:
            database = get_database(int(raw_database_id))
            required_db_type = str(BLAST_PROGRAMS[program]["db_type"])
            if database.db_type != required_db_type:
                raise ValueError(
                    f"{database.display_name} is {database.db_type}, but {program} requires {required_db_type}."
                )
            if database.status != "available":
                raise ValueError(f"{database.display_name} is marked as {database.status}.")

            result = run_blast(
                sequence=sequence,
                database=database.db_prefix_path,
                program=program,
                output_format=output_format,
                sensitivity_preset=sensitivity_preset,
                task=None if apoe_probe_preset else request.form.get("task") or None,
                evalue=request.form.get("evalue") or None,
                max_target_seqs=request.form.get("max_target_seqs") or None,
                word_size=request.form.get("word_size") or None,
                perc_identity="100" if apoe_probe_preset else request.form.get("perc_identity") or None,
                timeout_seconds=request.form.get("timeout_seconds") or None,
            )
            hits = filter_apoe_exact_probe_hits(result.hits) if apoe_probe_preset else result.hits
            saved_result = replace(result, hits=hits)
            run_id = save_result(saved_result)
            total_runtime_seconds += result.runtime_seconds
            total_hits += len(hits)
            query_count = result.query_count
            query_total_length = result.query_total_length
            database_results.append(
                {
                    "database_id": database.id,
                    "display_name": database.display_name,
                    "db_prefix_path": database.db_prefix_path,
                    "database_total_bytes": result.database_total_bytes,
                    "database_size_label": format_bytes(result.database_total_bytes),
                    "returncode": result.returncode,
                    "runtime_seconds": result.runtime_seconds,
                    "estimated_runtime_seconds": result.estimated_runtime_seconds,
                    "estimated_runtime_low_seconds": result.estimated_runtime_low_seconds,
                    "estimated_runtime_high_seconds": result.estimated_runtime_high_seconds,
                    "estimated_runtime_note": result.estimated_runtime_note,
                    "hit_count": len(hits),
                    "hits": hits,
                    "run_id": run_id,
                    "error": "",
                }
            )
        except Exception as exc:
            display_name = f"Database {raw_database_id}"
            try:
                display_name = get_database(int(raw_database_id)).display_name
            except Exception:
                pass
            database_results.append(
                {
                    "database_id": raw_database_id,
                    "display_name": display_name,
                    "db_prefix_path": "",
                    "database_total_bytes": 0,
                    "database_size_label": "unknown",
                    "returncode": "",
                    "runtime_seconds": "",
                    "estimated_runtime_seconds": None,
                    "estimated_runtime_low_seconds": None,
                    "estimated_runtime_high_seconds": None,
                    "estimated_runtime_note": "",
                    "hit_count": 0,
                    "hits": [],
                    "run_id": "",
                    "error": str(exc),
                }
            )

    payload = {
        "program": program,
        "output_format": output_format,
        "sensitivity_preset": sensitivity_preset,
        "query_count": query_count,
        "query_total_length": query_total_length,
        "total_runtime_seconds": total_runtime_seconds,
        "total_hits": total_hits,
        "apoe_probe_preset": apoe_probe_preset,
        "hit_filter": APOE_EXACT_MATCH_FILTER if apoe_probe_preset else "",
        "database_results": database_results,
    }
    batch_id = save_batch_result(payload)
    payload["batch_id"] = batch_id
    return render_template(
        "batch_results.html",
        batch=payload,
        error=None,
        format_duration=format_duration,
    )


@app.get("/databases")
def databases_page():
    """Render the database registry and database-management forms."""
    try:
        ensure_demo_databases()
        databases = list_databases()
        registry_error = None
    except Exception as exc:
        databases = []
        registry_error = str(exc)

    return render_template(
        "databases.html",
        categories=sorted(DB_CATEGORIES),
        db_types=sorted(DB_TYPES),
        databases=databases,
        error=request.args.get("error") or registry_error,
        message=request.args.get("message", ""),
    )


@app.get("/sra")
def sra_page():
    """Render local SRA discovery and pilot-database controls."""
    try:
        projects = discover_sra_projects()
        sra_blast_prefix_count = sum(len(project.blast_prefixes) for project in projects)
        error = request.args.get("error", "")
    except Exception as exc:
        projects = []
        sra_blast_prefix_count = 0
        error = request.args.get("error") or str(exc)

    toolkit_bin = sra_toolkit_bin()
    return render_template(
        "sra.html",
        projects=projects,
        sra_blast_prefix_count=sra_blast_prefix_count,
        sra_roots=[str(path) for path in configured_sra_roots()],
        sra_toolkit_bin=str(toolkit_bin) if toolkit_bin else "",
        message=request.args.get("message", ""),
        error=error,
    )


@app.post("/sra/register-db")
def register_sra_database_route():
    """Register an existing SRA-derived BLAST database prefix."""
    try:
        database = register_sra_blast_database(
            accession=request.form.get("accession", ""),
            db_prefix_path=request.form.get("db_prefix_path", ""),
        )
    except Exception as exc:
        return redirect_to_sra(error=str(exc))
    return redirect_to_sra(message=f"Registered {database.display_name}.")


@app.post("/sra/register-all-db")
def register_all_sra_databases_route():
    """Register every discovered SRA-derived BLAST database prefix."""
    registered = 0
    errors = []
    try:
        projects = discover_sra_projects()
    except Exception as exc:
        return redirect_to_sra(error=str(exc))

    for project in projects:
        for prefix in project.blast_prefixes:
            try:
                register_sra_blast_database(
                    accession=project.accession,
                    db_prefix_path=prefix,
                )
                registered += 1
            except Exception as exc:
                errors.append(f"{project.accession}: {exc}")

    if errors:
        return redirect_to_sra(
            message=f"Registered or updated {registered} SRA BLAST database prefix(es).",
            error="; ".join(errors[:5]),
        )
    if registered == 0:
        return redirect_to_sra(message="No discovered SRA BLAST databases were available to register.")
    return redirect_to_sra(message=f"Registered or updated {registered} SRA BLAST database prefix(es).")


@app.post("/sra/register-selected-db")
def register_selected_sra_databases_route():
    """Register the SRA-derived BLAST database prefixes selected in the workbench."""
    selected_databases = request.form.getlist("selected_db")
    if not selected_databases:
        return redirect_to_sra(error="Select at least one discovered SRA BLAST database to register.")

    registered = 0
    errors = []
    for selected_database in selected_databases:
        try:
            accession, prefix = selected_database.split("||", 1)
        except ValueError:
            errors.append("Skipped a selected database with an invalid form value.")
            continue

        try:
            register_sra_blast_database(
                accession=accession,
                db_prefix_path=prefix,
            )
            registered += 1
        except Exception as exc:
            errors.append(f"{accession}: {exc}")

    if errors:
        return redirect_to_sra(
            message=f"Registered or updated {registered} selected SRA BLAST database prefix(es).",
            error="; ".join(errors[:5]),
        )
    return redirect_to_sra(
        message=f"Registered or updated {registered} selected SRA BLAST database prefix(es)."
    )


@app.post("/sra/create-pilot")
def create_sra_pilot_route():
    """Create a small sampled BLAST database from an existing SRA FASTA file."""
    try:
        max_records = int(request.form.get("max_records", "1000"))
        database = create_pilot_database_from_fasta(
            accession=request.form.get("accession", ""),
            source_fasta_path=request.form.get("source_fasta_path", ""),
            max_records=max_records,
        )
    except Exception as exc:
        return redirect_to_sra(error=str(exc))
    return redirect_to_sra(message=f"Created pilot database {database.display_name}.")


@app.post("/sra/convert-pilot")
def convert_sra_pilot_route():
    """Convert a limited number of spots from local SRA to pilot FASTA."""
    try:
        max_spots = int(request.form.get("max_spots", "1000"))
        fasta_path = convert_sra_to_pilot_fasta(
            accession=request.form.get("accession", ""),
            sra_path=request.form.get("sra_path", ""),
            max_spots=max_spots,
        )
    except Exception as exc:
        return redirect_to_sra(error=str(exc))
    return redirect_to_sra(message=f"Created pilot FASTA: {fasta_path}")


@app.post("/databases/register")
def register_database_route():
    """Register an existing BLAST database prefix without modifying its files."""
    try:
        register_existing_database(
            display_name=request.form.get("display_name", ""),
            db_type=request.form.get("db_type", ""),
            db_prefix_path=request.form.get("db_prefix_path", ""),
            source_fasta_path=request.form.get("source_fasta_path") or None,
            description=request.form.get("description", ""),
            category=request.form.get("category", "custom"),
            notes=request.form.get("notes", ""),
        )
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(message="Database registered.")


@app.post("/databases/create")
def create_database_route():
    """Create BLAST index files from a FASTA file, then register the result."""
    try:
        create_database_from_fasta(
            display_name=request.form.get("display_name", ""),
            db_type=request.form.get("db_type", ""),
            source_fasta_path=request.form.get("source_fasta_path", ""),
            db_prefix_path=request.form.get("db_prefix_path") or None,
            description=request.form.get("description", ""),
            category=request.form.get("category", "custom"),
            notes=request.form.get("notes", ""),
        )
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(message="Database created with makeblastdb.")


@app.post("/databases/<int:database_id>/verify")
def verify_database_route(database_id: int):
    """Refresh one database's status by asking BLAST+ for its metadata."""
    try:
        database = verify_database(database_id)
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(
        message=f"{database.display_name} is marked as {database.status}."
    )


@app.post("/databases/verify-all")
def verify_all_databases_route():
    """Refresh every registered database status."""
    try:
        for database in list_databases():
            verify_database(database.id)
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(message="All registered databases were verified.")


@app.post("/databases/<int:database_id>/remove")
def remove_database_route(database_id: int):
    """Remove a database from the registry only; BLAST index files remain."""
    try:
        remove_database(database_id)
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(
        message="Database removed from the registry. BLAST files were not deleted."
    )


if __name__ == "__main__":
    app.run(host=FLASK_HOST, port=flask_port(), debug=False)
