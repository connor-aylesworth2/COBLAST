"""Flask entry point for the local COBLAST+ web interface.

This module keeps the web layer deliberately thin: routes collect form input,
delegate BLAST/database work to helper modules, and render templates with the
objects those helpers return.
"""

import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from time import perf_counter, time

from flask import (
    Flask,
    Response,
    abort,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)

from apoe_summary import apoe_probe_query_ids, build_apoe_probe_summary
from assembler import default_assembler
from contig_id import identify_contigs, reprobe_and_reassemble
from etol_summary import (
    ETOL_NET_FILTER,
    build_etol_matrix,
    build_etol_probe_summary,
    etol_control_query_ids,
    etol_preset_form_field,
    etol_preset_is_microbial,
    etol_preset_keys,
    etol_preset_label,
    etol_preset_options,
    etol_preset_records,
    etol_search_fasta,
    etol_search_query_ids,
    group_read_ids_by_taxon,
    sort_results_by_condition,
)
from etol_validation import compute_confusion
from design_matrix import DesignMatrixError, parse_design_matrix
from human_filter import extract_reads, filter_human_hits
from blast_runner import (
    BLAST_PROGRAMS,
    run_blast,
    run_blast_probe_panel,
    run_jobs_concurrently,
    validate_fasta_input,
)
from database_registry import (
    DB_CATEGORIES,
    DB_TYPES,
    create_database_from_fasta,
    ensure_demo_databases,
    get_database,
    list_databases,
    register_existing_database,
    remove_database,
    remove_missing_databases,
    verify_database,
)
from result_store import (
    apoe_summary_rows_as_delimited,
    batch_rows_as_delimited,
    batch_summary_rows_as_delimited,
    etol_confusion_rows_as_delimited,
    etol_contigs_as_fasta,
    etol_matrix_payload,
    etol_probe_counts_as_delimited,
    etol_summary_rows_as_delimited,
    load_batch_result,
    load_result,
    result_rows_as_delimited,
    save_batch_result,
    save_result,
)
from config import (
    FLASK_HOST,
    allocate_batch_resources,
    default_thread_count,
    flask_port,
    resource_path,
)
from database_size import database_storage_bytes, format_bytes
from sra_workflow import (
    configured_sra_roots,
    convert_sra_to_pilot_fasta,
    create_pilot_database_from_fasta,
    discover_sra_projects,
    register_sra_blast_database,
    source_fasta_for_blast_prefix,
    sra_toolkit_bin,
)


app = Flask(
    __name__,
    # resource_path works both from source and from a PyInstaller bundle.
    template_folder=str(resource_path("templates")),
    static_folder=str(resource_path("static")),
)

# Live progress for the synchronous batch route. The batch POST blocks while it
# runs, but Flask's dev server is threaded, so a separate poll request can read
# how many databases have finished. One local user and short-lived entries, so a
# plain dict under a lock is enough.
# ponytail: in-memory dict; needs a real store only if this ever serves >1 user.
_batch_progress: dict[str, dict] = {}
_batch_progress_lock = threading.Lock()


def _set_batch_stage(job_id: str, label: str, stage: str) -> None:
    """Record that database ``label`` has entered pipeline ``stage`` (live status).

    ``since`` lets the waiting page show how long a database has sat in one stage,
    which is the whole point for diagnosing slow steps (e.g. single-threaded CAP3).
    """
    if not job_id:
        return
    with _batch_progress_lock:
        entry = _batch_progress.get(job_id)
        if entry is not None:
            entry["stages"][label] = {"stage": stage, "since": time()}


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
APOE_PROBE_QUERY_IDS = apoe_probe_query_ids()


def numeric_hit_value(hit: dict[str, str], key: str) -> float | None:
    """Parse a numeric hit field that may be blank or formatted text."""
    try:
        return float(hit.get(key, ""))
    except (TypeError, ValueError):
        return None


def filter_exact_probe_hits(
    hits: list[dict[str, str]], probe_query_ids: set[str]
) -> list[dict[str, str]]:
    """Keep only exact probe matches (100% identity and coverage) after parsing.

    Used by the APOE genotyping preset, where only full-length exact matches are
    meaningful. The eToL presets use ``filter_net_probe_hits`` instead.
    """
    exact_hits = []
    for hit in hits:
        percent_identity = numeric_hit_value(hit, "pident")
        query_coverage = numeric_hit_value(hit, "qcovs")
        if (
            hit.get("qseqid", "") in probe_query_ids
            and percent_identity == 100.0
            and query_coverage == 100.0
        ):
            exact_hits.append(hit)
    return exact_hits


# Hu, Haas & Lathe 2022 gate the first-search net on E-value < 0.01 before
# counting (Abundance_ToL.py for the microbial probes, Abundance_count.py for the
# controls), so partial/mismatched rRNA matches are kept but statistically
# insignificant ones (down to BLAST's default e-value of 10) are not.
ETOL_EVALUE_THRESHOLD = 0.01


def filter_net_probe_hits(
    hits: list[dict[str, str]], probe_query_ids: set[str]
) -> list[dict[str, str]]:
    """Restrict eToL "net" hits to the active panel and gate them on E-value.

    Following Hu, Haas & Lathe 2022, the eToL panels cast a permissive net:
    default megablast with no identity or coverage filter, keeping partial and
    imperfect rRNA matches for the secondary human filter and cross-probe
    de-duplication to adjudicate. A hit is kept only if it belongs to the active
    panel and scores E-value < 0.01 (the paper's net cutoff, applied in
    Abundance_ToL.py); a hit whose E-value cannot be parsed is dropped.
    """
    kept = []
    for hit in hits:
        if hit.get("qseqid", "") not in probe_query_ids:
            continue
        evalue = numeric_hit_value(hit, "evalue")
        if evalue is None or evalue >= ETOL_EVALUE_THRESHOLD:
            continue
        kept.append(hit)
    return kept


def deduplicate_reads_to_best_probe(
    hits: list[dict[str, str]]
) -> tuple[list[dict[str, str]], int]:
    """Allocate each matched read to a single probe (Hu, Haas & Lathe 2022).

    A read (``sseqid``) recovered by more than one probe is counted only once,
    against the probe showing the highest sequence similarity. "Highest
    similarity" is ranked by bitscore, then percent identity, then query
    coverage, with the probe id as a final deterministic tie-break. Returns the
    de-duplicated hits (input order preserved) and the number of duplicate hits
    dropped.
    """
    best_by_read: dict[str, dict[str, str]] = {}
    for hit in hits:
        read_id = hit.get("sseqid", "")
        if not read_id:
            continue
        incumbent = best_by_read.get(read_id)
        if incumbent is None or _similarity_rank(hit) > _similarity_rank(incumbent):
            best_by_read[read_id] = hit

    winners = {id(hit) for hit in best_by_read.values()}
    kept = [hit for hit in hits if id(hit) in winners]
    return kept, len(hits) - len(kept)


def _similarity_rank(hit: dict[str, str]) -> tuple[float, float, float, str]:
    """Rank key for choosing a read's best probe: higher tuple wins.

    Ordered by bitscore, then percent identity, then query coverage. The probe id
    (``qseqid``) is the final element so the choice is a deterministic argmax
    independent of hit order; on a full score tie the lexically greater probe id
    wins (an arbitrary but stable rule).
    """
    return (
        numeric_hit_value(hit, "bitscore") or 0.0,
        numeric_hit_value(hit, "pident") or 0.0,
        numeric_hit_value(hit, "qcovs") or 0.0,
        hit.get("qseqid", ""),
    )


def count_control_reads(
    control_hits: list[dict[str, str]], control_query_ids: frozenset[str]
) -> dict[str, int]:
    """Tally host-normalization control reads, one read per best control probe.

    Mirrors Abundance_count.py section 2 (Hu, Haas & Lathe 2022): the control
    probe hits are de-duplicated so each read is allocated to the single control
    probe it matches best, then counted per probe (every control probe, including
    zeros). The E-value < 0.01 net gate is applied upstream by
    ``filter_net_probe_hits``. These counts are the host-cell normalization
    denominator, so counting one read against several redundant control probes
    would bias it.
    """
    deduped, _ = deduplicate_reads_to_best_probe(control_hits)
    counts = {probe_id: 0 for probe_id in control_query_ids}
    for hit in deduped:
        probe_id = hit.get("qseqid", "")
        if probe_id in counts:
            counts[probe_id] += 1
    return counts


def summarize_human_filter_warnings(database_results: list[dict]) -> str:
    """Roll per-database human-filter degradations into one batch-level warning.

    The human filter conservatively KEEPS any hit it cannot actually check --
    unrecoverable reads (no id-indexed DB and no readable source FASTA) or a
    failed human-genome BLAST -- and records why in each result's ``note``. Those
    notes only render in the per-database detail table, so a run that silently
    filtered nothing looks identical to a clean "0 removed" in the summary card
    (exactly how a mis-bundled build hides a broken filter). Surfacing the notes
    at batch level makes "couldn't run" impossible to mistake for "found none".
    """
    notes: list[str] = []
    for row in database_results:
        note = ((row.get("human_filter") or {}).get("note") or "").strip()
        if note and note not in notes:
            notes.append(note)
    return " ".join(notes)


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
        databases=databases,
        registry_error=registry_error,
        default_threads=default_thread_count(),
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
            task=request.form.get("task") or None,
            evalue=request.form.get("evalue") or None,
            max_target_seqs=request.form.get("max_target_seqs") or None,
            word_size=request.form.get("word_size") or None,
            perc_identity=request.form.get("perc_identity") or None,
            timeout_seconds=request.form.get("timeout_seconds") or None,
            num_threads=request.form.get("num_threads") or None,
        )
    except Exception as exc:
        return render_template(
            "results.html",
            error=str(exc),
            result=None,
            format_bytes=format_bytes,
        ), 400

    # Persist a JSON copy so the results page can offer CSV/TSV downloads.
    run_id = save_result(result)
    return render_template(
        "results.html",
        error=None,
        result=result,
        run_id=run_id,
        format_bytes=format_bytes,
    )


def _delimited_download(load, to_rows, filename_stem, file_format, guard=None):
    """Shared CSV/TSV download: validate format, load, optional guard, then stream.

    ``load`` is a zero-arg callable returning the saved payload (deferred so a bad
    extension 404s before any disk read); ``guard(payload)`` returns falsy to 404
    (e.g. the batch is not the expected preset); ``to_rows`` renders the payload
    as ``to_rows(payload, delimiter=...)``.
    """
    if file_format not in {"csv", "tsv"}:
        abort(404)
    try:
        payload = load()
    except FileNotFoundError:
        abort(404)
    if guard is not None and not guard(payload):
        abort(404)
    delimiter = "," if file_format == "csv" else "\t"
    mimetype = "text/csv" if file_format == "csv" else "text/tab-separated-values"
    return Response(
        to_rows(payload, delimiter=delimiter),
        mimetype=mimetype,
        headers={
            "Content-Disposition": f"attachment; filename={filename_stem}.{file_format}"
        },
    )


def _error_result_row(display_name, error, database_id=""):
    """Placeholder batch row for a database that failed before producing hits."""
    return {
        "database_id": database_id,
        "display_name": display_name,
        "db_prefix_path": "",
        "database_total_bytes": 0,
        "database_size_label": "unknown",
        "returncode": "",
        "runtime_seconds": "",
        "hit_count": 0,
        "hits": [],
        "run_id": "",
        "human_filter": None,
        "error": error,
    }


@app.get("/results/<run_id>.<file_format>")
def download_results(run_id: str, file_format: str):
    """Download a saved result table as CSV or TSV."""
    return _delimited_download(
        lambda: load_result(run_id),
        result_rows_as_delimited,
        f"blast_results_{run_id}",
        file_format,
    )


@app.get("/batch-results/<batch_id>.<file_format>")
def download_batch_results(batch_id: str, file_format: str):
    """Download saved batch results as CSV or TSV."""
    return _delimited_download(
        lambda: load_batch_result(batch_id),
        batch_rows_as_delimited,
        f"batch_blast_results_{batch_id}",
        file_format,
    )


@app.get("/batch-results/<batch_id>/summary.<file_format>")
def download_batch_summary(batch_id: str, file_format: str):
    """Download the Batch Summary panel statistics as CSV or TSV."""
    return _delimited_download(
        lambda: load_batch_result(batch_id),
        batch_summary_rows_as_delimited,
        f"batch_summary_{batch_id}",
        file_format,
    )


@app.get("/batch-results/<batch_id>/apoe-summary.<file_format>")
def download_apoe_summary(batch_id: str, file_format: str):
    """Download APOE probe count summaries as CSV or TSV."""
    return _delimited_download(
        lambda: load_batch_result(batch_id),
        apoe_summary_rows_as_delimited,
        f"apoe_probe_summary_{batch_id}",
        file_format,
        guard=lambda data: data.get("apoe_probe_preset"),
    )


@app.get("/batch-results/<batch_id>/etol-summary.<file_format>")
def download_etol_summary(batch_id: str, file_format: str):
    """Download eToL per-species exact-hit count summaries as CSV or TSV."""
    return _delimited_download(
        lambda: load_batch_result(batch_id),
        etol_summary_rows_as_delimited,
        f"etol_species_summary_{batch_id}",
        file_format,
        guard=lambda data: data.get("etol_probe_preset"),
    )


@app.get("/batch-results/<batch_id>/etol-probe-counts.<file_format>")
def download_etol_probe_counts(batch_id: str, file_format: str):
    """Download full per-probe exact-hit count data (every probe) as CSV or TSV."""
    return _delimited_download(
        lambda: load_batch_result(batch_id),
        etol_probe_counts_as_delimited,
        f"etol_probe_counts_{batch_id}",
        file_format,
        guard=lambda data: data.get("etol_probe_preset"),
    )


@app.get("/batch-results/<batch_id>/etol-confusion.<file_format>")
def download_etol_confusion(batch_id: str, file_format: str):
    """Download the eToL-V confusion matrix (per-cell) as CSV or TSV."""
    return _delimited_download(
        lambda: load_batch_result(batch_id),
        etol_confusion_rows_as_delimited,
        f"etol_v_confusion_{batch_id}",
        file_format,
        guard=lambda data: data.get("etol_preset_key") == "etol_v",
    )


@app.get("/batch-results/<batch_id>/etol-matrix.json")
def etol_matrix_json(batch_id: str):
    """Serve the plot-ready eToL hit matrix (rows x samples) for the heatmap."""
    try:
        batch_data = load_batch_result(batch_id)
    except FileNotFoundError:
        abort(404)

    if not batch_data.get("etol_probe_preset"):
        abort(404)

    level = request.args.get("level", "species")
    if level not in {"species", "probe"}:
        level = "species"
    return jsonify(etol_matrix_payload(batch_data, level=level))


@app.get("/batch-results/<batch_id>/etol-contigs.fasta")
def download_etol_contigs(batch_id: str):
    """Download all assembled eToL contigs for a batch as a multi-FASTA file."""
    try:
        batch_data = load_batch_result(batch_id)
    except FileNotFoundError:
        abort(404)

    if not batch_data.get("etol_probe_preset"):
        abort(404)

    body = etol_contigs_as_fasta(batch_data)
    filename = f"etol_contigs_{batch_id}.fasta"
    return Response(
        body,
        mimetype="text/plain",
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
        database_options=options,
        etol_preset_options=etol_preset_options(),
        error=request.args.get("error") or registry_error,
        message=request.args.get("message", ""),
    )


@app.get("/design-matrix-template.csv")
def design_matrix_template():
    """Download a starter design-matrix CSV in the strict ``sample,condition`` format."""
    body = (
        "sample,condition\n"
        "SRR21676099,AD\n"
        "SRR21676105,CONTROL\n"
        "SRR21676101,AD/LBD\n"
        "SRR21676126,AD/VaD\n"
    )
    return Response(
        body,
        mimetype="text/csv",
        headers={
            "Content-Disposition": "attachment; filename=design_matrix_template.csv"
        },
    )


@app.post("/batch-blast")
def run_batch_blast_route():
    """Run one query against multiple registered local databases concurrently."""
    sequence = request.form.get("sequence", "")
    uploaded_query = request.files.get("sequence_file")
    if uploaded_query and uploaded_query.filename:
        sequence = uploaded_query.read().decode("utf-8-sig")

    program = request.form.get("program", "blastn")
    database_ids = request.form.getlist("database_ids")
    apoe_probe_preset = request.form.get("apoe_probe_preset") == "1"
    # Exactly one exact-match probe preset can run at a time. The UI enforces
    # this, but the server resolves it deterministically too: pick the first
    # selected eToL preset, and an eToL preset takes precedence over APOE so the
    # stored query sets are never concatenated.
    etol_preset_key = next(
        (key for key in etol_preset_keys() if request.form.get(etol_preset_form_field(key)) == "1"),
        None,
    )
    if etol_preset_key:
        apoe_probe_preset = False
    exact_probe_preset = apoe_probe_preset or etol_preset_key is not None

    # Optional design matrix: an uploaded CSV/TSV that explicitly maps each
    # sample to a condition label for the eToL heatmap, overriding the unreliable
    # name-guessing regex. Parsed and validated up front so a malformed file is
    # rejected before the long BLAST run rather than after it. Only meaningful for
    # the eToL presets (the panels that render the heatmap).
    design_matrix_index: dict | None = None
    design_matrix_upload = request.files.get("design_matrix_file")
    if etol_preset_key and design_matrix_upload and design_matrix_upload.filename:
        try:
            design_matrix_index = parse_design_matrix(
                design_matrix_upload.read().decode("utf-8-sig"),
                filename=design_matrix_upload.filename,
            )
        except DesignMatrixError as exc:
            return (
                render_template(
                    "batch.html",
                    blast_programs=BLAST_PROGRAMS,
                    database_options=database_options(),
                    etol_preset_options=etol_preset_options(),
                    error=f"Design matrix error: {exc}",
                    message="",
                ),
                400,
            )

    # probe_query_ids is the full set of query ids searched (used to restrict
    # hits to the panel). For the microbial eToL presets the search also includes
    # the housekeeping control probes, so control_query_ids splits those out for
    # separate host-cell normalization; microbial_query_ids is the remainder.
    probe_query_ids: set[str] = set()
    control_query_ids: frozenset[str] = frozenset()
    if apoe_probe_preset:
        sequence = APOE_PROBE_FASTA
        probe_query_ids = APOE_PROBE_QUERY_IDS
    elif etol_preset_key:
        sequence = etol_search_fasta(etol_preset_key)
        probe_query_ids = set(etol_search_query_ids(etol_preset_key))
        if etol_preset_is_microbial(etol_preset_key):
            control_query_ids = etol_control_query_ids(etol_preset_key)
    if exact_probe_preset:
        # Exact-match/net probe presets always run BLASTN so the preset hit
        # filters below can be applied consistently.
        program = "blastn"

    # Secondary human filter: drop matched patient reads that also hit the human
    # genome. Only meaningful for the microbial eToL presets (the APOE/eToL
    # Control panels are human by design, so filtering human reads is nonsensical).
    human_filter_requested = request.form.get("human_filter") == "1"
    human_filter_active = (
        human_filter_requested
        and etol_preset_key is not None
        and etol_preset_is_microbial(etol_preset_key)
    )
    human_db = None
    if human_filter_active:
        try:
            human_db = get_database(int(request.form.get("human_filter_db_id", "")))
        except Exception:
            human_db = None
        if human_db is None or human_db.db_type != "nucl" or human_db.status != "available":
            return render_template(
                "batch.html",
                blast_programs=BLAST_PROGRAMS,
                database_options=database_options(),
                etol_preset_options=etol_preset_options(),
                error="Select an available nucleotide human-genome database for the secondary human filter.",
                message="",
            ), 400

    # Contig assembly (eToL re-probing input): assemble the reads each species'
    # probes recovered into longer contigs for species-level identification, as
    # Hu, Haas & Lathe 2022 do with CAP3. Only meaningful for the microbial eToL
    # presets. The CAP3 binary is optional, so a requested-but-unavailable
    # assembler degrades to a clear note rather than failing the run.
    assemble_contigs_requested = request.form.get("assemble_contigs") == "1"
    assemble_contigs_active = (
        assemble_contigs_requested
        and etol_preset_key is not None
        and etol_preset_is_microbial(etol_preset_key)
    )
    assembler = default_assembler()
    contig_assembly_available = assemble_contigs_active and assembler.is_available()
    contig_assembly_unavailable = assemble_contigs_active and not contig_assembly_available

    # Contig species identification (Hu, Haas & Lathe 2022): once contigs are
    # assembled, BLAST each against a local reference rRNA DB for its closest
    # homolog (species call) and re-BLAST each taxon's contigs against that
    # taxon's own reads to count near-100%-identity reads (confirmed abundance).
    # Requires contig assembly to be on plus an available nucleotide reference DB.
    # The paper's third step -- re-probing the library with key contigs -- is not
    # implemented yet.
    identify_contigs_requested = request.form.get("identify_contigs") == "1"
    identify_contigs_active = identify_contigs_requested and assemble_contigs_active
    reference_db = None
    if identify_contigs_active:
        try:
            reference_db = get_database(int(request.form.get("species_id_db_id", "")))
        except Exception:
            reference_db = None
        if reference_db is None or reference_db.db_type != "nucl" or reference_db.status != "available":
            return render_template(
                "batch.html",
                blast_programs=BLAST_PROGRAMS,
                database_options=database_options(),
                etol_preset_options=etol_preset_options(),
                error="Select an available nucleotide reference database (e.g. a SILVA rRNA database) for contig species identification.",
                message="",
            ), 400

    # Contig re-probing (Hu, Haas & Lathe 2022, Box 3): one round of using each
    # taxon's top contigs as probes against the SAME patient library to pull more
    # reads, then re-assembling. Searches the patient DB already in the batch (no
    # extra DB needed) and runs between assembly and identification, so naming and
    # confirmed abundance reflect the extended contigs. Requires assembly to be on.
    reprobe_requested = request.form.get("reprobe_contigs") == "1"
    reprobe_active = reprobe_requested and assemble_contigs_active

    if not database_ids:
        return render_template(
            "batch.html",
            blast_programs=BLAST_PROGRAMS,
            database_options=database_options(),
            etol_preset_options=etol_preset_options(),
            error="Choose at least one database for the batch run.",
            message="",
        ), 400

    # Pull every request-bound value out before spawning workers: Flask's
    # `request` is bound to this thread's context and must not be touched from
    # the concurrent worker threads below.
    task_value = None if exact_probe_preset else (request.form.get("task") or None)
    evalue_value = request.form.get("evalue") or None
    max_target_seqs_value = None if exact_probe_preset else (request.form.get("max_target_seqs") or None)
    word_size_value = request.form.get("word_size") or None
    perc_identity_value = None if exact_probe_preset else (request.form.get("perc_identity") or None)
    timeout_value = request.form.get("timeout_seconds") or None
    # Client-generated id so the waiting page can poll this batch's progress.
    job_id = request.form.get("job_id", "")

    # Validate the shared query once up front rather than inside every worker.
    # A bad query then fails fast with a single message instead of one identical
    # error per database, and the workers reuse the parsed result.
    try:
        prevalidated_query = validate_fasta_input(
            sequence, expected_type=str(BLAST_PROGRAMS[program]["query_type"])
        )
    except Exception as exc:
        return render_template(
            "batch.html",
            blast_programs=BLAST_PROGRAMS,
            database_options=database_options(),
            etol_preset_options=etol_preset_options(),
            error=str(exc),
            message="",
        ), 400

    # Benchmarks showed concurrency across patient databases scales far better
    # than -num_threads within one search, so split the core budget into
    # concurrent workers (most of it) plus a small per-job thread count.
    workers, threads_per_job = allocate_batch_resources(
        len(database_ids), requested_workers=request.form.get("batch_workers") or None
    )

    def run_single_database(raw_database_id):
        """Run one database's search and shape its result row (never raises)."""
        try:
            database = get_database(int(raw_database_id))
            required_db_type = str(BLAST_PROGRAMS[program]["db_type"])
            if database.db_type != required_db_type:
                raise ValueError(
                    f"{database.display_name} is {database.db_type}, but {program} requires {required_db_type}."
                )
            if database.status != "available":
                raise ValueError(f"{database.display_name} is marked as {database.status}.")

            _set_batch_stage(
                job_id,
                database.display_name,
                "Net BLAST search" if etol_preset_key else "BLAST search",
            )
            if etol_preset_key:
                # eToL panels run megablast over the whole panel (fast on
                # whole-SRA databases); a probe that cannot seed megablast is
                # silently dropped.
                result = run_blast_probe_panel(
                    panel_fasta=sequence,
                    database=database.db_prefix_path,
                    timeout_seconds=timeout_value,
                    num_threads=threads_per_job,
                )
            else:
                result = run_blast(
                    sequence=sequence,
                    database=database.db_prefix_path,
                    program=program,
                    # APOE exact-match and regular batch runs: run_blast applies
                    # the exact-match overrides (blastn-short, 100% identity and
                    # coverage, uncapped targets) when exact_match_probe is set.
                    task=task_value,
                    evalue=evalue_value,
                    max_target_seqs=max_target_seqs_value,
                    word_size=word_size_value,
                    perc_identity=perc_identity_value,
                    timeout_seconds=timeout_value,
                    num_threads=threads_per_job,
                    exact_match_probe=exact_probe_preset,
                    prevalidated_query=prevalidated_query,
                )
            # APOE genotyping keeps only full-length exact matches; the eToL
            # panels keep the paper's permissive net (default megablast, no
            # identity or coverage filter). Non-preset batch runs keep every hit.
            control_counts: dict[str, int] = {}
            if apoe_probe_preset:
                hits = filter_exact_probe_hits(result.hits, probe_query_ids)
            elif etol_preset_key:
                panel_hits = filter_net_probe_hits(result.hits, probe_query_ids)
                # Split the housekeeping control hits out BEFORE any human
                # filtering. Control reads (PGK1/hNSE) are human by design and
                # must never be human-filtered, or the host-cell normalization
                # denominator would be destroyed. They are de-duplicated to their
                # best control probe and counted per probe (including zeros) for
                # normalization, not as microbial species.
                control_hits = [
                    hit for hit in panel_hits
                    if hit.get("qseqid", "") in control_query_ids
                ]
                hits = [
                    hit for hit in panel_hits
                    if hit.get("qseqid", "") not in control_query_ids
                ]
                control_counts = count_control_reads(control_hits, control_query_ids)
            else:
                hits = result.hits
            # Secondary human filter runs on the microbial net hits only; a
            # failure here keeps the eToL hits unfiltered rather than discarding
            # the run.
            # Per-phase timing so a slow post-BLAST step (usually the human
            # filter's read BLAST or CAP3) shows up in the console after one run,
            # instead of bisecting by toggling features across 35-minute runs.
            phase_seconds = {"human_filter": 0.0, "dedup": 0.0, "assembly": 0.0, "reprobe": 0.0, "identification": 0.0}
            human_filter_stats = None
            if human_filter_active and human_db is not None:
                _set_batch_stage(job_id, database.display_name, "Secondary human filter")
                _phase_start = perf_counter()
                try:
                    hits, human_filter_stats = filter_human_hits(
                        hits,
                        db_prefix_path=database.db_prefix_path,
                        source_fasta_path=database.source_fasta_path,
                        human_db_prefix_path=human_db.db_prefix_path,
                        num_threads=threads_per_job,
                    )
                except Exception as exc:
                    human_filter_stats = {
                        "method": "error",
                        "reads_total": 0,
                        "reads_checked": 0,
                        "reads_unresolved": 0,
                        "human_reads": 0,
                        "hits_removed": 0,
                        "note": f"Human filter error: {exc}",
                    }
                phase_seconds["human_filter"] = perf_counter() - _phase_start
            # Cross-probe de-duplication (Hu, Haas & Lathe 2022): after human
            # removal, allocate each remaining read to the single probe with the
            # highest similarity so a read recovered by several redundant probes
            # is counted once. Runs last, as the paper specifies.
            dedup_removed = 0
            if etol_preset_key:
                _set_batch_stage(job_id, database.display_name, "De-duplicating reads")
                _phase_start = perf_counter()
                hits, dedup_removed = deduplicate_reads_to_best_probe(hits)
                phase_seconds["dedup"] = perf_counter() - _phase_start
            # Contig assembly (re-probing input): group the final non-human,
            # de-duplicated reads by species and assemble each group's reads into
            # contigs. The read sequences are recovered with the same helper the
            # human filter uses. A backend failure is recorded as a note rather
            # than dropping the database from the batch.
            contigs_by_species: dict[str, list[dict]] = {}
            contig_note = ""
            extract_seconds = 0.0
            cap3_seconds = 0.0
            # Kept in function scope so the identification step below can reuse
            # them; populated only when assembly actually runs.
            reads_by_taxon: dict[str, list[str]] = {}
            all_reads: dict[str, str] = {}
            if contig_assembly_available:
                _phase_start = perf_counter()
                try:
                    reads_by_taxon = group_read_ids_by_taxon(hits)
                    total_taxa = len(reads_by_taxon)
                    # Pull every species' reads out of the database in a SINGLE
                    # pass, then assemble from memory. The per-species version
                    # re-ran extract_reads (a blastdbcmd subprocess whose FASTA
                    # output is parsed in pure Python, plus a whole-file FASTA
                    # fallback) once per species, and that read extraction -- not
                    # CAP3 -- is what pinned a single core: its Python parsing
                    # holds the GIL, which the assembly thread pool cannot spread.
                    # One blastdbcmd call (one DB open, one parse) replaces N and
                    # leaves CAP3's subprocess as the only per-species work, which
                    # is what actually parallelizes across the shared pool.
                    all_read_ids = sorted(
                        {rid for ids in reads_by_taxon.values() for rid in ids}
                    )
                    _set_batch_stage(
                        job_id,
                        database.display_name,
                        f"Extracting {len(all_read_ids)} reads for assembly",
                    )
                    _extract_start = perf_counter()
                    all_reads, _method = extract_reads(
                        database.db_prefix_path,
                        database.source_fasta_path,
                        all_read_ids,
                    )
                    extract_seconds = perf_counter() - _extract_start

                    assembled = 0
                    progress_lock = threading.Lock()

                    def assemble_taxon(taxon, read_ids):
                        nonlocal assembled
                        reads = {
                            rid: all_reads[rid] for rid in read_ids if rid in all_reads
                        }
                        contigs = assembler.assemble(reads)
                        with progress_lock:
                            assembled += 1
                            _set_batch_stage(
                                job_id,
                                database.display_name,
                                f"Assembling contigs (CAP3): {assembled}/{total_taxa}",
                            )
                        return taxon, [contig.to_dict() for contig in contigs]

                    _cap3_start = perf_counter()
                    futures = [
                        assembly_pool.submit(assemble_taxon, taxon, read_ids)
                        for taxon, read_ids in reads_by_taxon.items()
                    ]
                    for future in futures:
                        # A failed taxon's exception surfaces here; let it propagate
                        # so the whole step degrades to a note, as the serial loop
                        # did, rather than silently dropping that database's contigs.
                        taxon, contig_dicts = future.result()
                        if contig_dicts:
                            contigs_by_species[taxon] = contig_dicts
                    cap3_seconds = perf_counter() - _cap3_start
                except Exception as exc:
                    contig_note = f"Contig assembly error: {exc}"
                phase_seconds["assembly"] = perf_counter() - _phase_start
            # Contig re-probing (one round): use each taxon's top contigs as
            # probes against this patient DB to pull more reads, then re-assemble.
            # Runs before identification so naming/abundance use the longer
            # contigs. Reuses the run's human DB when the human filter is on. A
            # failure degrades to a note, like assembly.
            reprobe_stats: dict[str, int] = {}
            reprobe_note = ""
            if reprobe_active and contig_assembly_available and contigs_by_species:
                _set_batch_stage(job_id, database.display_name, "Re-probing with key contigs")
                _phase_start = perf_counter()
                try:
                    reprobe_stats = reprobe_and_reassemble(
                        contigs_by_species,
                        reads_by_taxon,
                        all_reads,
                        patient_db_prefix=database.db_prefix_path,
                        source_fasta_path=database.source_fasta_path,
                        assembler=assembler,
                        num_threads=threads_per_job,
                        human_db_prefix=(
                            human_db.db_prefix_path
                            if (human_filter_active and human_db is not None)
                            else None
                        ),
                        assembly_pool=assembly_pool,
                    )
                except Exception as exc:
                    reprobe_note = f"Re-probing error: {exc}"
                phase_seconds["reprobe"] = perf_counter() - _phase_start
            # Contig species identification + confirmed abundance: BLAST the
            # assembled contigs against the reference rRNA DB (closest homolog)
            # and against each taxon's own reads (>= identity). Annotates the
            # contig dicts in place. A failure degrades to a note, like assembly.
            contig_identification: dict[str, dict] = {}
            contig_id_note = ""
            if identify_contigs_active and reference_db is not None and contigs_by_species:
                _set_batch_stage(
                    job_id,
                    database.display_name,
                    "Identifying contigs (species ID + confirmed abundance)",
                )
                _phase_start = perf_counter()
                try:
                    contig_identification = identify_contigs(
                        contigs_by_species,
                        reads_by_taxon,
                        all_reads,
                        reference_db_prefix=reference_db.db_prefix_path,
                        num_threads=threads_per_job,
                    )
                except Exception as exc:
                    contig_id_note = f"Contig identification error: {exc}"
                phase_seconds["identification"] = perf_counter() - _phase_start
            print(
                f"[etol-timing] {database.display_name}: "
                f"blast={result.runtime_seconds:.1f}s "
                f"human_filter={phase_seconds['human_filter']:.1f}s "
                f"dedup={phase_seconds['dedup']:.1f}s "
                f"assembly={phase_seconds['assembly']:.1f}s "
                f"reprobe={phase_seconds['reprobe']:.1f}s "
                f"identification={phase_seconds['identification']:.1f}s "
                f"(extract={extract_seconds:.1f}s cap3={cap3_seconds:.1f}s)",
                flush=True,
            )
            _set_batch_stage(job_id, database.display_name, "Saving results")
            run_id = save_result(replace(result, hits=hits))
            _set_batch_stage(job_id, database.display_name, "Done")
            row = {
                "database_id": database.id,
                "display_name": database.display_name,
                "db_prefix_path": database.db_prefix_path,
                "database_total_bytes": result.database_total_bytes,
                "database_size_label": format_bytes(result.database_total_bytes),
                "returncode": result.returncode,
                "runtime_seconds": result.runtime_seconds,
                "hit_count": len(hits),
                "hits": hits,
                "run_id": run_id,
                "human_filter": human_filter_stats,
                "etol_control_counts": control_counts,
                "etol_dedup_removed": dedup_removed,
                "contigs": contigs_by_species,
                "contig_count": sum(len(items) for items in contigs_by_species.values()),
                "contig_species_count": len(contigs_by_species),
                "contig_note": contig_note,
                "reprobe_new_reads": reprobe_stats.get("new_reads", 0),
                "reprobe_taxa": reprobe_stats.get("reprobed_taxa", 0),
                "reprobe_human_removed": reprobe_stats.get("human_removed", 0),
                "reprobe_note": reprobe_note,
                "contig_identification": contig_identification,
                "contig_id_note": contig_id_note,
                "contigs_identified": sum(
                    1
                    for items in contigs_by_species.values()
                    for contig in items
                    if contig.get("closest_homolog")
                ),
                "error": "",
            }
            return {
                "row": row,
                "runtime": result.runtime_seconds,
                "hits": len(hits),
                "query_count": result.query_count,
                "query_total_length": result.query_total_length,
            }
        except Exception as exc:
            display_name = f"Database {raw_database_id}"
            try:
                display_name = get_database(int(raw_database_id)).display_name
            except Exception:
                pass
            return {
                "row": _error_result_row(display_name, str(exc), database_id=raw_database_id),
                "runtime": 0.0,
                "hits": 0,
                "query_count": 0,
                "query_total_length": 0,
            }
        finally:
            # One database finished (success or error); advance the live counter
            # the waiting page polls.
            if job_id:
                with _batch_progress_lock:
                    entry = _batch_progress.get(job_id)
                    if entry:
                        entry["done"] += 1

    # Publish the total up front so the first poll shows a determinate bar.
    if job_id:
        with _batch_progress_lock:
            _batch_progress[job_id] = {"done": 0, "total": len(database_ids), "stages": {}}

    # Each BLAST search is a separate process, so threads give real parallelism.
    # One assembly pool is shared across every database for the whole batch and
    # capped at the core budget, so CAP3 work load-balances across databases:
    # a single-DB job uses every core for its species, and the last database
    # still assembling in a many-DB job claims the freed cores instead of
    # crawling through its species one core at a time. run_single_database
    # (a closure) submits its per-species assemblies here.
    # ponytail: assembly is capped at budget independently of the BLAST fan, so
    # the two can briefly sum above budget while some databases still search;
    # unify them under one scheduler only if that overlap measurably thrashes.
    wall_start = perf_counter()
    try:
        with ThreadPoolExecutor(max_workers=default_thread_count()) as assembly_pool:
            outcomes = run_jobs_concurrently(
                run_single_database,
                [{"raw_database_id": raw_database_id} for raw_database_id in database_ids],
                max_workers=workers,
            )
    finally:
        if job_id:
            with _batch_progress_lock:
                _batch_progress.pop(job_id, None)
    wall_clock_seconds = perf_counter() - wall_start

    database_results = []
    total_runtime_seconds = 0.0
    total_hits = 0
    query_count = 0
    query_total_length = 0
    for outcome in outcomes:
        if isinstance(outcome, Exception):
            # Defensive: run_single_database catches its own errors, but never
            # let a leaked exception drop a database from the results table.
            database_results.append(_error_result_row("Unknown database", str(outcome)))
            continue
        database_results.append(outcome["row"])
        total_runtime_seconds += outcome["runtime"]
        total_hits += outcome["hits"]
        if outcome["query_count"]:
            query_count = outcome["query_count"]
            query_total_length = outcome["query_total_length"]

    # Group samples by design-matrix condition once, before summaries/exports are
    # built and the batch is saved, so the heatmap, tables, and downloads share
    # one condition-ordered column layout.
    if design_matrix_index:
        database_results = sort_results_by_condition(database_results, design_matrix_index)

    hit_filter = ""
    if apoe_probe_preset:
        hit_filter = APOE_EXACT_MATCH_FILTER
    elif etol_preset_key:
        hit_filter = ETOL_NET_FILTER

    payload = {
        "program": program,
        "query_count": query_count,
        "query_total_length": query_total_length,
        "total_runtime_seconds": total_runtime_seconds,
        "wall_clock_seconds": wall_clock_seconds,
        "batch_workers": workers,
        "total_hits": total_hits,
        "apoe_probe_preset": apoe_probe_preset,
        "etol_probe_preset": etol_preset_key is not None,
        "etol_preset_key": etol_preset_key,
        "etol_preset_label": etol_preset_label(etol_preset_key) if etol_preset_key else "",
        "hit_filter": hit_filter,
        "human_filter_enabled": human_filter_active,
        "human_filter_db": human_db.display_name if human_db else "",
        "human_filter_hits_removed": sum(
            (result_row.get("human_filter") or {}).get("hits_removed", 0)
            for result_row in database_results
        ),
        "human_filter_warning": summarize_human_filter_warnings(database_results),
        "etol_normalized": etol_preset_key is not None
        and etol_preset_is_microbial(etol_preset_key),
        "etol_dedup_removed": sum(
            result_row.get("etol_dedup_removed", 0) for result_row in database_results
        ),
        "assemble_contigs": assemble_contigs_active,
        "contig_assembly_unavailable": contig_assembly_unavailable,
        "contig_count": sum(
            result_row.get("contig_count", 0) for result_row in database_results
        ),
        "identify_contigs": identify_contigs_active,
        "species_id_db": reference_db.display_name if reference_db else "",
        "contigs_identified": sum(
            result_row.get("contigs_identified", 0) for result_row in database_results
        ),
        "reprobe_contigs": reprobe_active,
        "reprobe_new_reads": sum(
            result_row.get("reprobe_new_reads", 0) for result_row in database_results
        ),
        "database_results": database_results,
    }
    # Persist the uploaded design matrix (if any) so the heatmap JSON endpoint --
    # which reloads the saved batch from disk -- can apply it for condition labels.
    if design_matrix_index:
        payload["design_matrix"] = design_matrix_index
    if apoe_probe_preset:
        payload["apoe_probe_summary"] = build_apoe_probe_summary(database_results)
    if etol_preset_key:
        payload["etol_probe_summary"] = build_etol_probe_summary(
            database_results, etol_preset_records(etol_preset_key)
        )
    # eToL-V runs get a confusion matrix vs the bundled eToL WGS ground truth
    # (the dissertation's Figure 9 fidelity check). Guarded so a validation
    # failure never takes down the results page.
    if etol_preset_key == "etol_v":
        try:
            matrix = build_etol_matrix(
                database_results,
                etol_preset_records("etol_v"),
                condition_index=design_matrix_index,
            )
            payload["etol_confusion"] = compute_confusion(matrix)
        except Exception as exc:  # pragma: no cover - defensive
            payload["etol_confusion"] = {"error": str(exc)}
    batch_id = save_batch_result(payload)
    payload["batch_id"] = batch_id
    return render_template(
        "batch_results.html",
        batch=payload,
        error=None,
    )


@app.get("/batch-progress/<job_id>")
def batch_progress(job_id: str):
    """Report how many databases of an in-flight batch have finished.

    Polled by the waiting page to fill its progress bar. Unknown ids (not started
    yet, or already cleaned up) report zeros so the client stays neutral.
    """
    with _batch_progress_lock:
        entry = _batch_progress.get(job_id)
        if not entry:
            return {"done": 0, "total": 0, "stages": []}
        return {
            "done": entry["done"],
            "total": entry["total"],
            "stages": [
                {"label": label, "stage": info["stage"], "since": info["since"]}
                for label, info in entry.get("stages", {}).items()
            ],
        }


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
        missing_database_count=sum(
            database.status == "missing" for database in databases
        ),
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
            source_fasta_path=request.form.get("source_fasta_path") or None,
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
                    source_fasta_path=source_fasta_for_blast_prefix(
                        prefix, project.fasta_files
                    )
                    or None,
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
    try:
        projects = discover_sra_projects()
    except Exception as exc:
        return redirect_to_sra(error=str(exc))
    source_fasta_by_prefix = {
        prefix: source_fasta_for_blast_prefix(prefix, project.fasta_files)
        for project in projects
        for prefix in project.blast_prefixes
    }
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
                source_fasta_path=source_fasta_by_prefix.get(prefix) or None,
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


@app.post("/databases/remove-missing")
def remove_missing_databases_route():
    """Remove every database currently marked missing from the registry."""
    try:
        removed_databases = remove_missing_databases()
    except Exception as exc:
        return redirect_to_databases(error=str(exc))

    removed_count = len(removed_databases)
    if removed_count == 0:
        return redirect_to_databases(message="No missing databases were found.")
    noun = "database" if removed_count == 1 else "databases"
    return redirect_to_databases(
        message=(
            f"Removed {removed_count} missing {noun} from the registry. "
            "BLAST files were not deleted."
        )
    )


def _parse_selected_database_ids() -> tuple[list[int], list[str]]:
    """Parse the ``selected_db`` checkboxes into integer ids, collecting errors."""
    database_ids: list[int] = []
    errors: list[str] = []
    for raw_id in request.form.getlist("selected_db"):
        try:
            database_ids.append(int(raw_id))
        except (TypeError, ValueError):
            errors.append(f"Skipped a selection with an invalid id: {raw_id!r}.")
    return database_ids, errors


@app.post("/databases/verify-selected")
def verify_selected_databases_route():
    """Refresh the status of each database selected on the registry page."""
    database_ids, errors = _parse_selected_database_ids()
    if not database_ids:
        return redirect_to_databases(error="Select at least one database to verify.")

    verified = 0
    for database_id in database_ids:
        try:
            verify_database(database_id)
            verified += 1
        except Exception as exc:
            errors.append(f"Database {database_id}: {exc}")

    noun = "database" if verified == 1 else "databases"
    message = f"Verified {verified} selected {noun}."
    if errors:
        return redirect_to_databases(message=message, error="; ".join(errors[:5]))
    return redirect_to_databases(message=message)


@app.post("/databases/remove-selected")
def remove_selected_databases_route():
    """Remove each selected database from the registry; BLAST files remain."""
    database_ids, errors = _parse_selected_database_ids()
    if not database_ids:
        return redirect_to_databases(error="Select at least one database to remove.")

    removed = 0
    for database_id in database_ids:
        try:
            remove_database(database_id)
            removed += 1
        except Exception as exc:
            errors.append(f"Database {database_id}: {exc}")

    noun = "database" if removed == 1 else "databases"
    message = (
        f"Removed {removed} selected {noun} from the registry. "
        "BLAST files were not deleted."
    )
    if errors:
        return redirect_to_databases(message=message, error="; ".join(errors[:5]))
    return redirect_to_databases(message=message)


@app.post("/databases/<int:database_id>/remove")
def remove_database_route(database_id: int):
    """Remove a database from the registry only; BLAST index files remain."""
    try:
        database = remove_database(database_id)
    except Exception as exc:
        return redirect_to_databases(error=str(exc))
    return redirect_to_databases(
        message=(
            f"{database.display_name} was removed from the registry. "
            "BLAST files were not deleted."
        )
    )


if __name__ == "__main__":
    app.run(host=FLASK_HOST, port=flask_port(), debug=False)
