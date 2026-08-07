"""Checks for the downstream figure script's parsing and scoring helpers.

Only the logic that can be wrong *quietly* is covered. `homolog_genus` earns a
test because its first implementation split the title on whitespace, returned
the accession, and scored every named genus in the published AD shortlist as
absent -- a wrong F9 that was indistinguishable from a real null result.

Fixtures are real `Closest homolog (contig)` strings from an EBB eToL Full run.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from plot_thesis import bh, cohens_d, homolog_genus, student_p  # noqa: E402


def test_homolog_genus_reads_the_second_to_last_rank():
    cases = {
        # Shortlist genera that the accession-splitting bug missed entirely.
        "CP023705.2380511.2382005 Bacteria;Pseudomonadota;Alphaproteobacteria;"
        "Sphingomonadales;Sphingomonadaceae;Sphingomonas;Sphingomonas melonis": "Sphingomonas",
        "JQ698928.1.1530 Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;"
        "Dikarya;Ascomycota;Saccharomycotina;Saccharomycetes;Saccharomycetales;"
        "Phaffomycetaceae;Komagataella;Komagataella pastoris": "Komagataella",
        # Last rank is not a binomial: splitting it on whitespace would give
        # "Sphingomonadaceae", the family, not the genus.
        "DQ490372.1.1450 Bacteria;Pseudomonadota;Alphaproteobacteria;"
        "Sphingomonadales;Sphingomonadaceae;Sphingomonas;"
        "Sphingomonadaceae bacterium KVD-unk-19": "Sphingomonas",
        "FPLS01033064.15.1481 Bacteria;Pseudomonadota;Alphaproteobacteria;"
        "Sphingomonadales;Sphingomonadaceae;Sphingomonas;metagenome": "Sphingomonas",
        # A land plant, which is what the Chloroplastida group actually carries.
        "LRBV01001970.46707.48505 Eukaryota;Archaeplastida;Chloroplastida;Charophyta;"
        "Phragmoplastophyta;Streptophyta;Embryophyta;Tracheophyta;Spermatophyta;"
        "Magnoliophyta;Fagales;Quercus;Quercus lobata": "Quercus",
    }
    for title, genus in cases.items():
        assert homolog_genus(title) == genus, title[:60]

    # Never returns the accession, whatever the title looks like.
    assert not homolog_genus("").startswith("LRBV")
    assert homolog_genus("no-lineage-here") == ""


def test_ported_statistics_match_closed_forms():
    # Two-sided t critical value at alpha = 0.05, df = 10.
    assert abs(student_p(2.228, 10) - 0.05) < 1e-4
    # The canonical Benjamini-Hochberg worked example; output stays monotone.
    adjusted = bh([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205])
    assert abs(adjusted[0] - 0.008) < 1e-6
    assert adjusted == sorted(adjusted)
    # Means 6 vs 2, pooled SD 1.
    assert abs(cohens_d([5, 6, 7], [1, 2, 3]) - 4.0) < 1e-9
