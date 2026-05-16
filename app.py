from pathlib import Path
import tempfile

from flask import Flask, render_template, request

from blast_runner import BLAST_PROGRAMS, run_blast


app = Flask(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_TOY_DB = Path(tempfile.gettempdir()) / "blast_flask_demo" / "db" / "toy_nt"


@app.get("/")
def index():
    return render_template(
        "index.html",
        blast_programs=BLAST_PROGRAMS,
        default_db=DEFAULT_TOY_DB,
    )


@app.post("/run-blast")
def run_blast_route():
    sequence = request.form.get("sequence", "")
    database = request.form.get("database", str(DEFAULT_TOY_DB))
    program = request.form.get("program", "blastn")
    output_format = request.form.get("output_format", "tabular")

    try:
        result = run_blast(
            sequence=sequence,
            database=database,
            program=program,
            output_format=output_format,
        )
    except Exception as exc:
        return render_template("results.html", error=str(exc), result=None), 400

    return render_template("results.html", error=None, result=result)


if __name__ == "__main__":
    app.run(debug=True)
