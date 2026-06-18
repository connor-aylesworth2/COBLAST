"""Pluggable contig assembly for the eToL re-probing workflow.

The eToL net (Hu, Haas & Lathe, *BMC Microbiology* 2022;22:317) recovers many
short rRNA reads per probe. Because each 64-mer probe detects a *cluster* of
related species rather than a unique one, the paper pins the exact species by
assembling a probe/species' matched reads into longer contigs (it used CAP3 and
EGassembler), then re-identifying and re-probing with those contigs.

This module provides that assembly step behind a small interface so the eToL
batch route depends on the abstraction "reads -> contigs" rather than on any one
assembler. :class:`Cap3Assembler` is the default backend (the tool the paper
used); the :class:`Assembler` protocol lets another engine -- or a no-op in
tests -- be substituted without touching the pipeline.

CAP3 is filename-driven: it reads one input FASTA and writes its results to
sibling files (``<input>.cap.contigs``, ``.cap.singlets``, ``.cap.ace`` ...), so
each assembly runs in its own temporary directory and concurrent per-species
assemblies (the batch route fans databases across a thread pool) cannot collide
on those output names.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable
import re
import subprocess
import tempfile

from config import cap3_exe


@dataclass(frozen=True)
class Contig:
    """One assembled contig plus how many reads were merged into it.

    ``num_reads`` is the read support reported by the assembler (0 when it could
    not be parsed); the re-probing workflow ranks contigs by it to find the
    "most abundant" sequences, as the paper does.
    """

    id: str
    sequence: str
    num_reads: int = 0

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable dict (batch results are stored as JSON)."""
        return {
            "id": self.id,
            "sequence": self.sequence,
            "num_reads": self.num_reads,
            "length": self.length,
        }


@runtime_checkable
class Assembler(Protocol):
    """Turns a set of reads into assembled contigs.

    Implementations must be safe to call from worker threads and must not raise
    for the ordinary "nothing assembled" outcome -- they return an empty list
    instead. A genuine backend failure (e.g. the assembler process erroring) may
    raise, so the caller can record a note and keep the rest of the run alive.
    """

    name: str

    def is_available(self) -> bool:
        """Return True when the backend can actually run (e.g. binary present)."""
        ...

    def assemble(self, reads: dict[str, str]) -> list[Contig]:
        """Assemble ``{read_id: sequence}`` into contigs (possibly empty)."""
        ...


def _reads_to_fasta(reads: dict[str, str]) -> str:
    """Render ``{read_id: sequence}`` as FASTA text for the assembler input."""
    return "".join(f">{read_id}\n{sequence}\n" for read_id, sequence in reads.items())


# CAP3's ACE ``CO`` (contig) lines look like ``CO <name> <nbases> <nreads> ...``;
# the fourth field is the read count, which we map back onto each contig.
_ACE_CONTIG_LINE = re.compile(r"^CO\s+(\S+)\s+\d+\s+(\d+)\b")


def _ace_read_counts(ace_text: str) -> dict[str, int]:
    """Map contig id -> number of reads, parsed from CAP3's ``.cap.ace`` output."""
    counts: dict[str, int] = {}
    for line in ace_text.splitlines():
        match = _ACE_CONTIG_LINE.match(line.strip())
        if match:
            counts[match.group(1)] = int(match.group(2))
    return counts


def parse_cap3_contigs(contigs_text: str, read_counts: dict[str, int] | None = None) -> list[Contig]:
    """Parse CAP3's ``.cap.contigs`` FASTA into :class:`Contig` objects.

    ``read_counts`` (from :func:`_ace_read_counts`) supplies per-contig read
    support; contigs missing from it default to 0.
    """
    read_counts = read_counts or {}
    contigs: list[Contig] = []
    name: str | None = None
    seq_parts: list[str] = []
    for line in contigs_text.splitlines():
        stripped = line.strip()
        if stripped.startswith(">"):
            if name is not None and seq_parts:
                contigs.append(Contig(name, "".join(seq_parts), read_counts.get(name, 0)))
            header = stripped[1:].strip()
            name = header.split()[0] if header else None
            seq_parts = []
        elif stripped:
            seq_parts.append(stripped)
    if name is not None and seq_parts:
        contigs.append(Contig(name, "".join(seq_parts), read_counts.get(name, 0)))
    return contigs


class Cap3Assembler:
    """CAP3 overlap-layout-consensus backend (the assembler used by the paper).

    ``overlap_length`` and ``overlap_identity_pct`` map to CAP3's ``-o`` and
    ``-p`` overlap cutoffs; left as ``None`` they use CAP3's own defaults so
    results match the paper's web-server runs unless deliberately tuned.
    """

    name = "CAP3"

    # Minimum reads worth assembling: a single read is, by definition, its own
    # singlet and produces no contig, so CAP3 is not invoked below this.
    MIN_READS = 2

    def __init__(
        self,
        overlap_length: str | None = None,
        overlap_identity_pct: str | None = None,
    ) -> None:
        self.overlap_length = overlap_length
        self.overlap_identity_pct = overlap_identity_pct

    def is_available(self) -> bool:
        """Return True when a CAP3 binary can be resolved."""
        try:
            cap3_exe()
        except FileNotFoundError:
            return False
        return True

    def assemble(self, reads: dict[str, str]) -> list[Contig]:
        """Assemble reads with CAP3, returning the contigs (empty if none form).

        Raises ``RuntimeError`` if CAP3 exits non-zero (a real failure, distinct
        from the ordinary outcome of every read staying a singlet, which returns
        an empty list).
        """
        usable = {read_id: sequence for read_id, sequence in reads.items() if sequence}
        if len(usable) < self.MIN_READS:
            return []

        exe = cap3_exe()  # FileNotFoundError here is a caller bug; gate on is_available().
        with tempfile.TemporaryDirectory(prefix="cap3_") as tmpdir:
            reads_path = Path(tmpdir) / "reads.fasta"
            reads_path.write_text(_reads_to_fasta(usable), encoding="utf-8")

            # CAP3 writes outputs next to (and named after) the input file, so we
            # run with cwd=tmpdir and reference the input by name.
            command = [str(exe), reads_path.name]
            if self.overlap_length:
                command += ["-o", str(self.overlap_length)]
            if self.overlap_identity_pct:
                command += ["-p", str(self.overlap_identity_pct)]

            completed = subprocess.run(
                command,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                check=False,
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    "CAP3 assembly failed: "
                    + ((completed.stderr or completed.stdout).strip() or "unknown error")
                )

            contigs_path = Path(tmpdir) / "reads.fasta.cap.contigs"
            if not contigs_path.exists():
                # Return code 0 but no contigs file: every read stayed a singlet.
                return []

            ace_path = Path(tmpdir) / "reads.fasta.cap.ace"
            read_counts = (
                _ace_read_counts(ace_path.read_text(encoding="utf-8", errors="replace"))
                if ace_path.exists()
                else {}
            )
            return parse_cap3_contigs(
                contigs_path.read_text(encoding="utf-8", errors="replace"),
                read_counts,
            )


def default_assembler() -> Assembler:
    """Return the default contig-assembly backend (CAP3, CAP3 defaults)."""
    return Cap3Assembler()
