from pathlib import Path

from flask import Flask, redirect, render_template, request, url_for

from blast_runner import BLAST_PROGRAMS, run_blast
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


app = Flask(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent


def redirect_to_databases(message: str = "", error: str = ""):
    params = {}
    if message:
        params["message"] = message
    if error:
        params["error"] = error
    return redirect(url_for("databases_page", **params))


@app.get("/")
def index():
    try:
        ensure_demo_databases()
        databases = list_databases()
        registry_error = None
    except Exception as exc:
        databases = []
        registry_error = str(exc)

    return render_template(
        "index.html",
        blast_programs=BLAST_PROGRAMS,
        databases=databases,
        registry_error=registry_error,
    )


@app.post("/run-blast")
def run_blast_route():
    sequence = request.form.get("sequence", "")
    program = request.form.get("program", "blastn")
    database_id = request.form.get("database_id", "")
    output_format = request.form.get("output_format", "tabular")

    try:
        if program not in BLAST_PROGRAMS:
            raise ValueError(f"Unsupported BLAST program: {program}")
        if not database_id:
            raise ValueError("Choose a registered BLAST database.")

        database = get_database(int(database_id))
        required_db_type = str(BLAST_PROGRAMS[program]["db_type"])
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

        result = run_blast(
            sequence=sequence,
            database=database.db_prefix_path,
            program=program,
            output_format=output_format,
        )
    except Exception as exc:
        return render_template("results.html", error=str(exc), result=None), 400

    return render_template("results.html", error=None, result=result)


@app.get("/databases")
def databases_page():
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
    try:
        database = verify_database(database_id)
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(
        message=f"{database.display_name} is marked as {database.status}."
    )


@app.post("/databases/verify-all")
def verify_all_databases_route():
    try:
        for database in list_databases():
            verify_database(database.id)
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(message="All registered databases were verified.")


@app.post("/databases/<int:database_id>/remove")
def remove_database_route(database_id: int):
    try:
        remove_database(database_id)
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(
        message="Database removed from the registry. BLAST files were not deleted."
    )


if __name__ == "__main__":
    app.run(debug=True)
