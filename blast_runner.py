"""Validation, command construction, and parsing for local BLAST+ searches.

The Flask routes pass plain form values into this module. The functions here
normalize FASTA input, enforce safe local-only options, call the appropriate
BLAST+ executable, and convert stdout into table rows for the interface.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import subprocess
import tempfile
from time import perf_counter
from typing import Any

from config import DISALLOWED_BLAST_OPTIONS, REMOTE_BLAST_ENABLED, blast_exe
from runtime_estimator import database_storage_bytes, estimate_blast_runtime_seconds

try:
    from Bio import SearchIO, SeqIO
except ImportError as exc:  # pragma: no cover - exercised only when dependency is absent
    SearchIO = None
    SeqIO = None
    BIOPYTHON_IMPORT_ERROR = exc
else:
    BIOPYTHON_IMPORT_ERROR = None


OUTFMT6_FIELDS = [
    "qseqid",
    "sseqid",
    "stitle",
    "pident",
    "length",
    "qcovs",
    "evalue",
    "bitscore",
]
ALLOWED_BLASTN_TASKS = {"blastn", "blastn-short", "dc-megablast", "megablast"}

# Each program defines the query type users must submit and the database type
# that must be selected. The UI reads this same mapping to filter databases.
BLAST_PROGRAMS = {
    "blastn": {
        "label": "BLASTN",
        "description": "nucleotide query vs nucleotide database",
        "query_type": "nucleotide",
        "db_type": "nucl",
        "default_task": "blastn-short",
        "allowed_tasks": ALLOWED_BLASTN_TASKS,
    },
    "blastp": {
        "label": "BLASTP",
        "description": "protein query vs protein database",
        "query_type": "protein",
        "db_type": "prot",
        "default_task": None,
        "allowed_tasks": set(),
    },
    "blastx": {
        "label": "BLASTX",
        "description": "translated nucleotide query vs protein database",
        "query_type": "nucleotide",
        "db_type": "prot",
        "default_task": None,
        "allowed_tasks": set(),
    },
    "tblastn": {
        "label": "TBLASTN",
        "description": "protein query vs translated nucleotide database",
        "query_type": "protein",
        "db_type": "nucl",
        "default_task": None,
        "allowed_tasks": set(),
    },
}
BLAST_OUTPUT_FORMATS = {
    # Format 6 is tabular; naming the columns keeps SearchIO parsing predictable.
    "tabular": "6 " + " ".join(OUTFMT6_FIELDS),
    "xml": "5",
}
FAST_TIMEOUT_SECONDS = 300
DEFAULT_TIMEOUT_SECONDS = 600
SENSITIVE_TIMEOUT_SECONDS = 900
SENSITIVITY_PRESETS = {
    # Presets keep common searches one-click while advanced fields can override
    # individual values later in build_blast_parameters.
    "standard": {
        "label": "Standard",
        "description": "Balanced default for routine sequence checks.",
        "evalue": "10",
        "max_target_seqs": "50",
        "timeout_seconds": DEFAULT_TIMEOUT_SECONDS,
    },
    "sensitive": {
        "label": "Sensitive",
        "description": "Keeps weaker candidate matches for review.",
        "evalue": "100",
        "max_target_seqs": "100",
        "timeout_seconds": SENSITIVE_TIMEOUT_SECONDS,
    },
    "fast": {
        "label": "Fast",
        "description": "Returns a smaller hit list for quick checks.",
        "evalue": "10",
        "max_target_seqs": "10",
        "timeout_seconds": FAST_TIMEOUT_SECONDS,
    },
}
NUCLEOTIDE_ALPHABET = set("ACGTRYSWKMBDHVNU")
PROTEIN_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")
MAX_FASTA_RECORDS = 100
MAX_TOTAL_SEQUENCE_LENGTH = 5_000_000
FASTA_LINE_WIDTH = 80
MAX_TARGET_SEQS_LIMIT = 10_000
TIMEOUT_SECONDS_LIMIT = 3_600


@dataclass(frozen=True)
class FastaRecordSummary:
    """Small per-record summary used in the results page."""

    id: str
    length: int


@dataclass(frozen=True)
class FastaValidationResult:
    """Normalized FASTA text plus metadata gathered during validation."""

    fasta: str
    sequence_type: str
    records: list[FastaRecordSummary]
    total_length: int


@dataclass(frozen=True)
class BlastResult:
    """Complete outcome of one BLAST run, including command and parsed hits."""

    returncode: int
    hits: list[dict[str, str]]
    stdout: str
    stderr: str
    command: list[str]
    database_path: str
    database_total_bytes: int
    output_format: str
    program: str
    runtime_seconds: float
    estimated_runtime_seconds: float | None
    estimated_runtime_low_seconds: float | None
    estimated_runtime_high_seconds: float | None
    estimated_runtime_note: str
    query_type: str
    query_count: int
    query_total_length: int
    sensitivity_preset: str
    parameters: dict[str, str]


def require_biopython() -> None:
    """Fail early with an installation hint if Biopython is unavailable."""
    if SearchIO is None or SeqIO is None:
        raise RuntimeError(
            "Biopython is required for FASTA validation and BLAST result parsing. "
            "Install dependencies with: python -m pip install -r requirements.txt"
        ) from BIOPYTHON_IMPORT_ERROR


def wrap_sequence(sequence: str, width: int = FASTA_LINE_WIDTH) -> str:
    """Wrap sequence text at a conventional FASTA line width."""
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def coerce_to_fasta_text(sequence: str) -> str:
    """Accept either FASTA text or a bare sequence and return FASTA text."""
    cleaned = sequence.strip()
    if not cleaned:
        raise ValueError("Enter a FASTA sequence.")

    normalized = cleaned.replace("\r\n", "\n").replace("\r", "\n")
    if normalized.lstrip().startswith(">"):
        return normalized

    # A pasted sequence without a header is still valid input for the interface.
    raw_sequence = "".join(normalized.split())
    if not raw_sequence:
        raise ValueError("Enter a sequence with at least one residue or base.")
    return f">query\n{raw_sequence}"


def normalize_fasta_lines(fasta_text: str) -> str:
    """Clean spacing, require headers, and uppercase sequence lines."""
    normalized_lines: list[str] = []
    seen_header = False

    for line_number, line in enumerate(fasta_text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            header = stripped[1:].strip()
            if not header:
                raise ValueError(f"FASTA header on line {line_number} is empty.")
            normalized_lines.append(f">{header}")
            seen_header = True
            continue
        if not seen_header:
            raise ValueError(
                f"Sequence data appears before the first FASTA header on line {line_number}."
            )
        normalized_lines.append("".join(stripped.split()).upper())

    if not normalized_lines:
        raise ValueError("Enter a FASTA sequence.")

    return "\n".join(normalized_lines) + "\n"


def invalid_characters(sequence: str, expected_type: str) -> set[str]:
    """Return residues/bases that do not belong to the selected query type."""
    if expected_type == "nucleotide":
        return set(sequence) - NUCLEOTIDE_ALPHABET
    if expected_type == "protein":
        return set(sequence) - PROTEIN_ALPHABET
    raise ValueError(f"Unsupported query sequence type: {expected_type}")


def validate_fasta_input(
    sequence: str,
    expected_type: str = "nucleotide",
) -> FastaValidationResult:
    """Parse and validate FASTA before BLAST sees it.

    BLAST accepts a broad range of input, but the interface benefits from clear
    errors and predictable normalized text. This also catches accidentally using
    protein characters with nucleotide programs, or vice versa.
    """
    require_biopython()
    if expected_type not in {"nucleotide", "protein"}:
        raise ValueError(f"Unsupported query sequence type: {expected_type}")

    fasta_text = normalize_fasta_lines(coerce_to_fasta_text(sequence))
    records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
    if not records:
        raise ValueError("No FASTA records could be parsed from the query.")
    if len(records) > MAX_FASTA_RECORDS:
        raise ValueError(
            f"Too many FASTA records ({len(records)}). "
            f"The current prototype accepts up to {MAX_FASTA_RECORDS} records per run."
        )

    seen_ids: set[str] = set()
    summaries: list[FastaRecordSummary] = []
    normalized_records: list[str] = []
    total_length = 0

    for record_number, record in enumerate(records, start=1):
        # IDs need to be stable because they become the qseqid shown in tables.
        record_id = record.id.strip()
        if not record_id or record_id == "<unknown id>":
            raise ValueError(f"FASTA record {record_number} does not have a usable ID.")
        if record_id in seen_ids:
            raise ValueError(f"Duplicate FASTA record ID: {record_id}")
        seen_ids.add(record_id)

        seq = str(record.seq).replace(" ", "").replace("\t", "").upper()
        if not seq:
            raise ValueError(f"FASTA record {record_id} has no sequence.")
        if "-" in seq or "." in seq:
            raise ValueError(
                f"FASTA record {record_id} contains gap characters. "
                "Remove '-' or '.' before running BLAST."
            )

        invalid = invalid_characters(seq, expected_type)
        if invalid:
            chars = ", ".join(sorted(invalid))
            raise ValueError(
                f"FASTA record {record_id} contains characters that are not valid "
                f"for a {expected_type} query: {chars}"
            )

        if expected_type == "nucleotide":
            # U is allowed for pasted RNA-like inputs, then normalized to DNA.
            seq = seq.replace("U", "T")

        total_length += len(seq)
        if total_length > MAX_TOTAL_SEQUENCE_LENGTH:
            raise ValueError(
                f"Query contains {total_length:,} total bases/residues. "
                f"The current prototype limit is {MAX_TOTAL_SEQUENCE_LENGTH:,}."
            )

        description = record.description.strip() or record_id
        normalized_records.append(f">{description}\n{wrap_sequence(seq)}")
        summaries.append(FastaRecordSummary(id=record_id, length=len(seq)))

    return FastaValidationResult(
        fasta="\n".join(normalized_records) + "\n",
        sequence_type=expected_type,
        records=summaries,
        total_length=total_length,
    )


def validate_fasta(sequence: str, expected_type: str = "nucleotide") -> str:
    """Compatibility wrapper for callers that only need normalized FASTA."""
    return validate_fasta_input(sequence, expected_type=expected_type).fasta


def require_searchio() -> None:
    """Alias used by parsing functions to make their dependency explicit."""
    require_biopython()


def format_float(value: Any, decimals: int) -> str:
    """Render optional numeric values for table cells."""
    if value is None or value == "":
        return ""
    return f"{float(value):.{decimals}f}"


def format_evalue(value: Any) -> str:
    """Render e-values consistently for the results table."""
    if value is None:
        return ""
    number = float(value)
    if number == 0:
        return "0.0"
    return f"{number:.2e}"


def percent_identity(hsp: Any) -> float | None:
    """Read percent identity from SearchIO, with a manual fallback."""
    if getattr(hsp, "ident_pct", None) is not None:
        return float(hsp.ident_pct)

    ident_num = getattr(hsp, "ident_num", None)
    aln_span = getattr(hsp, "aln_span", None)
    if ident_num is None or not aln_span:
        return None
    return (float(ident_num) / float(aln_span)) * 100


def query_coverage(qresult: Any, hit: Any, hsp: Any) -> float | None:
    """Read or calculate query coverage as a percentage."""
    if getattr(hit, "query_coverage", None) is not None:
        return float(hit.query_coverage)

    query_span = getattr(hsp, "query_span", None)
    query_length = getattr(qresult, "seq_len", None)
    if query_span is None or not query_length:
        return None
    return (float(query_span) / float(query_length)) * 100


def subject_title(hit: Any, hsp: Any) -> str:
    """Choose the most descriptive subject label available from SearchIO."""
    for candidate in (
        getattr(hit, "title", None),
        getattr(hit, "description", None),
        getattr(hsp, "hit_description", None),
        getattr(hsp, "hit_id", None),
        getattr(hit, "id", None),
    ):
        if candidate and candidate != "<unknown description>":
            return str(candidate)
    return ""


def searchio_results_to_hits(qresults: Iterable[Any]) -> list[dict[str, str]]:
    """Flatten SearchIO query/hit/HSP objects into table-row dictionaries."""
    hits: list[dict[str, str]] = []
    for qresult in qresults:
        for hit in qresult:
            for hsp in hit:
                hits.append(
                    {
                        "qseqid": getattr(hsp, "query_id", None) or qresult.id,
                        "sseqid": getattr(hsp, "hit_id", None) or hit.id,
                        "stitle": subject_title(hit, hsp),
                        "pident": format_float(percent_identity(hsp), 3),
                        "length": str(getattr(hsp, "aln_span", "")),
                        "qcovs": format_float(query_coverage(qresult, hit, hsp), 1),
                        "evalue": format_evalue(getattr(hsp, "evalue", None)),
                        "bitscore": format_float(getattr(hsp, "bitscore", None), 1),
                    }
                )
    return hits


def parse_blast_tabular(stdout: str) -> list[dict[str, str]]:
    """Parse BLAST format-6 stdout into result rows."""
    if not stdout.strip():
        return []

    hits: list[dict[str, str]] = []
    expected_column_count = len(OUTFMT6_FIELDS)
    for line_number, raw_line in enumerate(stdout.splitlines(), start=1):
        if not raw_line.strip():
            continue
        values = raw_line.rstrip("\n").split("\t")
        if len(values) != expected_column_count:
            raise ValueError(
                "Could not parse BLAST tabular output on line "
                f"{line_number}: expected {expected_column_count} columns, found {len(values)}."
            )
        row = dict(zip(OUTFMT6_FIELDS, values, strict=True))
        hits.append(
            {
                "qseqid": row["qseqid"],
                "sseqid": row["sseqid"],
                "stitle": row["stitle"],
                "pident": format_float(row["pident"], 3),
                "length": row["length"],
                "qcovs": format_float(row["qcovs"], 1),
                "evalue": format_evalue(row["evalue"]),
                "bitscore": format_float(row["bitscore"], 1),
            }
        )
    return hits


def parse_blast_xml(stdout: str) -> list[dict[str, str]]:
    """Parse BLAST XML stdout into result rows."""
    require_searchio()
    if not stdout.strip():
        return []
    qresults = SearchIO.parse(StringIO(stdout), "blast-xml")
    return searchio_results_to_hits(qresults)


def parse_blast_output(stdout: str, output_format: str = "tabular") -> list[dict[str, str]]:
    """Dispatch to the parser that matches the selected output format."""
    if output_format == "tabular":
        return parse_blast_tabular(stdout)
    if output_format == "xml":
        return parse_blast_xml(stdout)
    raise ValueError(f"Unsupported BLAST output format: {output_format}")


def enforce_local_blast_only(command: list[str]) -> None:
    """Prevent accidental remote BLAST usage from this local-only prototype."""
    if REMOTE_BLAST_ENABLED:
        raise RuntimeError("Remote BLAST cannot be enabled for this local interface.")

    command_options = {part.lower() for part in command}
    blocked_options = {option.lower() for option in DISALLOWED_BLAST_OPTIONS}
    used_blocked_options = sorted(command_options & blocked_options)
    if used_blocked_options:
        raise RuntimeError(
            "Remote BLAST is disabled for this local interface. "
            f"Blocked option(s): {', '.join(used_blocked_options)}"
        )


def optional_text(value: str | None) -> str | None:
    """Normalize optional form fields so blank strings behave like missing data."""
    if value is None:
        return None
    cleaned = value.strip()
    return cleaned or None


def parse_positive_float(name: str, value: str | None) -> str | None:
    """Validate a positive numeric option and return its original string value."""
    cleaned = optional_text(value)
    if cleaned is None:
        return None
    try:
        number = float(cleaned)
    except ValueError as exc:
        raise ValueError(f"{name} must be a number.") from exc
    if number <= 0:
        raise ValueError(f"{name} must be greater than 0.")
    return cleaned


def parse_bounded_int(name: str, value: str | None, minimum: int, maximum: int) -> str | None:
    """Validate an integer option with inclusive minimum/maximum bounds."""
    cleaned = optional_text(value)
    if cleaned is None:
        return None
    try:
        number = int(cleaned)
    except ValueError as exc:
        raise ValueError(f"{name} must be a whole number.") from exc
    if number < minimum or number > maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}.")
    return str(number)


def preset_timeout_seconds(sensitivity_preset: str) -> str:
    """Return the timeout attached to the selected sensitivity preset."""
    if sensitivity_preset not in SENSITIVITY_PRESETS:
        allowed = ", ".join(SENSITIVITY_PRESETS)
        raise ValueError(f"Unsupported sensitivity preset: {sensitivity_preset}. Choose one of: {allowed}.")
    return str(SENSITIVITY_PRESETS[sensitivity_preset]["timeout_seconds"])


def parse_percent_identity(program: str, value: str | None) -> str | None:
    """Validate BLASTN-only percent identity filtering."""
    cleaned = optional_text(value)
    if cleaned is None:
        return None
    if program != "blastn":
        raise ValueError("Minimum percent identity is currently supported for BLASTN only.")
    try:
        number = float(cleaned)
    except ValueError as exc:
        raise ValueError("Minimum percent identity must be a number.") from exc
    if number < 0 or number > 100:
        raise ValueError("Minimum percent identity must be between 0 and 100.")
    return cleaned


def build_blast_parameters(
    *,
    program: str,
    sensitivity_preset: str,
    evalue: str | None,
    max_target_seqs: str | None,
    word_size: str | None,
    perc_identity: str | None,
) -> dict[str, str]:
    """Merge preset defaults with user-supplied advanced BLAST options."""
    if sensitivity_preset not in SENSITIVITY_PRESETS:
        allowed = ", ".join(SENSITIVITY_PRESETS)
        raise ValueError(f"Unsupported sensitivity preset: {sensitivity_preset}. Choose one of: {allowed}.")

    preset = SENSITIVITY_PRESETS[sensitivity_preset]
    # Start with a safe preset, then overlay validated advanced fields below.
    parameters = {
        "evalue": str(preset["evalue"]),
        "max_target_seqs": str(preset["max_target_seqs"]),
    }

    parsed_evalue = parse_positive_float("E-value", evalue)
    parsed_max_target_seqs = parse_bounded_int(
        "Maximum target sequences",
        max_target_seqs,
        1,
        MAX_TARGET_SEQS_LIMIT,
    )
    parsed_word_size = parse_bounded_int(
        "Word size",
        word_size,
        4 if program == "blastn" else 2,
        1_000,
    )
    parsed_perc_identity = parse_percent_identity(program, perc_identity)

    if parsed_evalue is not None:
        parameters["evalue"] = parsed_evalue
    if parsed_max_target_seqs is not None:
        parameters["max_target_seqs"] = parsed_max_target_seqs
    if parsed_word_size is not None:
        parameters["word_size"] = parsed_word_size
    if parsed_perc_identity is not None:
        parameters["perc_identity"] = parsed_perc_identity

    return parameters


def run_blast(
    sequence: str,
    database: str | Path,
    program: str = "blastn",
    timeout_seconds: int | str | None = None,
    task: str | None = None,
    output_format: str = "tabular",
    sensitivity_preset: str = "standard",
    evalue: str | None = None,
    max_target_seqs: str | None = None,
    word_size: str | None = None,
    perc_identity: str | None = None,
) -> BlastResult:
    """Run one local BLAST search and return both raw and parsed outputs."""
    if program not in BLAST_PROGRAMS:
        allowed = ", ".join(BLAST_PROGRAMS)
        raise ValueError(f"Unsupported BLAST program: {program}. Choose one of: {allowed}.")
    if output_format not in BLAST_OUTPUT_FORMATS:
        raise ValueError(f"Unsupported BLAST output format: {output_format}")
    timeout_value = str(timeout_seconds) if timeout_seconds is not None else None
    if optional_text(timeout_value) is None:
        timeout_value = preset_timeout_seconds(sensitivity_preset)
    # subprocess.run enforces this timeout, so keep it bounded for the UI.
    timeout = parse_bounded_int(
        "Timeout",
        timeout_value,
        1,
        TIMEOUT_SECONDS_LIMIT,
    )

    program_config = BLAST_PROGRAMS[program]
    default_task = program_config["default_task"]
    # Only BLASTN exposes task variants in this interface.
    selected_task = task if task is not None else default_task
    allowed_tasks = program_config["allowed_tasks"]
    if selected_task is not None and selected_task not in allowed_tasks:
        raise ValueError(f"Unsupported task for {program}: {selected_task}")
    parameters = build_blast_parameters(
        program=program,
        sensitivity_preset=sensitivity_preset,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        word_size=word_size,
        perc_identity=perc_identity,
    )

    query = validate_fasta_input(
        sequence,
        expected_type=str(program_config["query_type"]),
    )
    db_path = str(database)
    database_total_bytes = database_storage_bytes(db_path)
    runtime_estimate = estimate_blast_runtime_seconds(
        program=program,
        query_total_length=query.total_length,
        database_bytes=database_total_bytes,
        sensitivity_preset=sensitivity_preset,
    )

    with tempfile.TemporaryDirectory(prefix="blast_flask_") as tmpdir:
        query_path = Path(tmpdir) / "query.fasta"
        # BLAST+ expects a file path for -query, so the pasted/uploaded sequence
        # lives in a short-lived temporary FASTA file.
        query_path.write_text(query.fasta, encoding="utf-8")

        cmd = [
            str(blast_exe(program)),
            "-query",
            str(query_path),
            "-db",
            db_path,
            "-outfmt",
            BLAST_OUTPUT_FORMATS[output_format],
        ]
        if selected_task is not None:
            cmd.extend(["-task", selected_task])
        for parameter, value in parameters.items():
            cmd.extend([f"-{parameter}", value])
        enforce_local_blast_only(cmd)

        start = perf_counter()
        # Capture stdout/stderr so the interface can display parsed hits and
        # still expose BLAST diagnostics when a run returns warnings/errors.
        completed = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=int(timeout),
            check=False,
        )
        runtime_seconds = perf_counter() - start

    return BlastResult(
        returncode=completed.returncode,
        hits=parse_blast_output(completed.stdout, output_format)
        if completed.returncode == 0
        else [],
        stdout=completed.stdout,
        stderr=completed.stderr,
        command=cmd,
        database_path=db_path,
        database_total_bytes=database_total_bytes,
        output_format=output_format,
        program=program,
        runtime_seconds=runtime_seconds,
        estimated_runtime_seconds=runtime_estimate.seconds if runtime_estimate else None,
        estimated_runtime_low_seconds=runtime_estimate.low_seconds if runtime_estimate else None,
        estimated_runtime_high_seconds=runtime_estimate.high_seconds if runtime_estimate else None,
        estimated_runtime_note=runtime_estimate.note if runtime_estimate else "",
        query_type=query.sequence_type,
        query_count=len(query.records),
        query_total_length=query.total_length,
        sensitivity_preset=sensitivity_preset,
        parameters=parameters,
    )


def run_blastn(
    sequence: str,
    database: str | Path,
    timeout_seconds: int | str | None = None,
    task: str = "blastn-short",
    output_format: str = "tabular",
    sensitivity_preset: str = "standard",
) -> BlastResult:
    """Convenience wrapper retained for older blastn-only callers/tests."""
    return run_blast(
        sequence=sequence,
        database=database,
        program="blastn",
        timeout_seconds=timeout_seconds,
        task=task,
        output_format=output_format,
        sensitivity_preset=sensitivity_preset,
    )
