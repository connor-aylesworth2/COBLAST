"""Flask entry point for the local COBLAST+ web interface.

This module keeps the web layer deliberately thin: routes collect form input,
delegate BLAST/database work to helper modules, and render templates with the
objects those helpers return.
"""

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
from result_store import load_result, result_rows_as_delimited, save_result
from config import FLASK_HOST, flask_port, resource_path, resource_root


app = Flask(
    __name__,
    # resource_path works both from source and from a PyInstaller bundle.
    template_folder=str(resource_path("templates")),
    static_folder=str(resource_path("static")),
)

PROJECT_ROOT = resource_root()


def redirect_to_databases(message: str = "", error: str = ""):
    """Send users back to the database page with optional status text."""
    params = {}
    if message:
        params["message"] = message
    if error:
        params["error"] = error
    return redirect(url_for("databases_page", **params))


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
        return render_template("results.html", error=str(exc), result=None), 400

    # Persist a JSON copy so the results page can offer CSV/TSV downloads.
    run_id = save_result(result)
    return render_template("results.html", error=None, result=result, run_id=run_id)


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
