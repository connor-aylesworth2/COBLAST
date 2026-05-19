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
    "tabular": "6 " + " ".join(OUTFMT6_FIELDS),
    "xml": "5",
}
FAST_TIMEOUT_SECONDS = 300
DEFAULT_TIMEOUT_SECONDS = 600
SENSITIVE_TIMEOUT_SECONDS = 900
SENSITIVITY_PRESETS = {
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
    id: str
    length: int


@dataclass(frozen=True)
class FastaValidationResult:
    fasta: str
    sequence_type: str
    records: list[FastaRecordSummary]
    total_length: int


@dataclass(frozen=True)
class BlastResult:
    returncode: int
    hits: list[dict[str, str]]
    stdout: str
    stderr: str
    command: list[str]
    output_format: str
    program: str
    runtime_seconds: float
    query_type: str
    query_count: int
    query_total_length: int
    sensitivity_preset: str
    parameters: dict[str, str]


def require_biopython() -> None:
    if SearchIO is None or SeqIO is None:
        raise RuntimeError(
            "Biopython is required for FASTA validation and BLAST result parsing. "
            "Install dependencies with: python -m pip install -r requirements.txt"
        ) from BIOPYTHON_IMPORT_ERROR


def wrap_sequence(sequence: str, width: int = FASTA_LINE_WIDTH) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def coerce_to_fasta_text(sequence: str) -> str:
    cleaned = sequence.strip()
    if not cleaned:
        raise ValueError("Enter a FASTA sequence.")

    normalized = cleaned.replace("\r\n", "\n").replace("\r", "\n")
    if normalized.lstrip().startswith(">"):
        return normalized

    raw_sequence = "".join(normalized.split())
    if not raw_sequence:
        raise ValueError("Enter a sequence with at least one residue or base.")
    return f">query\n{raw_sequence}"


def normalize_fasta_lines(fasta_text: str) -> str:
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
    if expected_type == "nucleotide":
        return set(sequence) - NUCLEOTIDE_ALPHABET
    if expected_type == "protein":
        return set(sequence) - PROTEIN_ALPHABET
    raise ValueError(f"Unsupported query sequence type: {expected_type}")


def validate_fasta_input(
    sequence: str,
    expected_type: str = "nucleotide",
) -> FastaValidationResult:
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
    return validate_fasta_input(sequence, expected_type=expected_type).fasta


def require_searchio() -> None:
    require_biopython()


def format_float(value: Any, decimals: int) -> str:
    if value is None:
        return ""
    return f"{float(value):.{decimals}f}"


def format_evalue(value: Any) -> str:
    if value is None:
        return ""
    number = float(value)
    if number == 0:
        return "0.0"
    return f"{number:.2e}"


def percent_identity(hsp: Any) -> float | None:
    if getattr(hsp, "ident_pct", None) is not None:
        return float(hsp.ident_pct)

    ident_num = getattr(hsp, "ident_num", None)
    aln_span = getattr(hsp, "aln_span", None)
    if ident_num is None or not aln_span:
        return None
    return (float(ident_num) / float(aln_span)) * 100


def query_coverage(qresult: Any, hit: Any, hsp: Any) -> float | None:
    if getattr(hit, "query_coverage", None) is not None:
        return float(hit.query_coverage)

    query_span = getattr(hsp, "query_span", None)
    query_length = getattr(qresult, "seq_len", None)
    if query_span is None or not query_length:
        return None
    return (float(query_span) / float(query_length)) * 100


def subject_title(hit: Any, hsp: Any) -> str:
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
    require_searchio()
    if not stdout.strip():
        return []
    qresults = SearchIO.parse(StringIO(stdout), "blast-tab", fields=OUTFMT6_FIELDS)
    return searchio_results_to_hits(qresults)


def parse_blast_xml(stdout: str) -> list[dict[str, str]]:
    require_searchio()
    if not stdout.strip():
        return []
    qresults = SearchIO.parse(StringIO(stdout), "blast-xml")
    return searchio_results_to_hits(qresults)


def parse_blast_output(stdout: str, output_format: str = "tabular") -> list[dict[str, str]]:
    if output_format == "tabular":
        return parse_blast_tabular(stdout)
    if output_format == "xml":
        return parse_blast_xml(stdout)
    raise ValueError(f"Unsupported BLAST output format: {output_format}")


def enforce_local_blast_only(command: list[str]) -> None:
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
    if value is None:
        return None
    cleaned = value.strip()
    return cleaned or None


def parse_positive_float(name: str, value: str | None) -> str | None:
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
    if sensitivity_preset not in SENSITIVITY_PRESETS:
        allowed = ", ".join(SENSITIVITY_PRESETS)
        raise ValueError(f"Unsupported sensitivity preset: {sensitivity_preset}. Choose one of: {allowed}.")
    return str(SENSITIVITY_PRESETS[sensitivity_preset]["timeout_seconds"])


def parse_percent_identity(program: str, value: str | None) -> str | None:
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
    if sensitivity_preset not in SENSITIVITY_PRESETS:
        allowed = ", ".join(SENSITIVITY_PRESETS)
        raise ValueError(f"Unsupported sensitivity preset: {sensitivity_preset}. Choose one of: {allowed}.")

    preset = SENSITIVITY_PRESETS[sensitivity_preset]
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
    if program not in BLAST_PROGRAMS:
        allowed = ", ".join(BLAST_PROGRAMS)
        raise ValueError(f"Unsupported BLAST program: {program}. Choose one of: {allowed}.")
    if output_format not in BLAST_OUTPUT_FORMATS:
        raise ValueError(f"Unsupported BLAST output format: {output_format}")
    timeout_value = str(timeout_seconds) if timeout_seconds is not None else None
    if optional_text(timeout_value) is None:
        timeout_value = preset_timeout_seconds(sensitivity_preset)
    timeout = parse_bounded_int(
        "Timeout",
        timeout_value,
        1,
        TIMEOUT_SECONDS_LIMIT,
    )

    program_config = BLAST_PROGRAMS[program]
    default_task = program_config["default_task"]
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

    with tempfile.TemporaryDirectory(prefix="blast_flask_") as tmpdir:
        query_path = Path(tmpdir) / "query.fasta"
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
        output_format=output_format,
        program=program,
        runtime_seconds=runtime_seconds,
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
    return run_blast(
        sequence=sequence,
        database=database,
        program="blastn",
        timeout_seconds=timeout_seconds,
        task=task,
        output_format=output_format,
        sensitivity_preset=sensitivity_preset,
    )
