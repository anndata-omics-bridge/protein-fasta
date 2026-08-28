"""Tests for aggregate protein diagnostics."""

from protein_fasta.diagnostic_summary import summarize_protein_diagnostics
from protein_fasta.record import ProteinDiagnostics, ProteinRecord


def test_summary_preserves_overlaps_and_counts_changes() -> None:
    diagnostics = [
        ProteinDiagnostics(
            protein=ProteinRecord("REV_CON_P1", "first", "AC?"),
            raw_header="REV_CON_P1 first",
            identifier_namespace="unknown",
            classifications=frozenset({"decoy", "contaminant"}),
            upper_cased=True,
            stop_stripped=True,
            illegal_residues="?",
        ),
        ProteinDiagnostics(
            protein=ProteinRecord("P2", None, "M??"),
            raw_header="P2",
            identifier_namespace="uniprot_bare",
            classifications=frozenset(),
            upper_cased=False,
            stop_stripped=False,
            illegal_residues="??",
        ),
    ]

    summary = summarize_protein_diagnostics(diagnostics)

    assert summary.proteins.n_sequences == 2
    assert summary.proteins.total_residues == 6
    assert summary.namespace_counts == {"unknown": 1, "uniprot_bare": 1}
    assert summary.classification_counts == {"contaminant": 1, "decoy": 1}
    assert summary.classification_combination_counts == {
        (): 1,
        ("contaminant", "decoy"): 1,
    }
    assert summary.upper_cased_count == 1
    assert summary.stop_stripped_count == 1
    assert summary.illegal_sequence_count == 2
    assert summary.illegal_residue_counts == {"?": 3}


def test_empty_summary_has_stable_empty_values() -> None:
    summary = summarize_protein_diagnostics([])

    assert summary.proteins.n_sequences == 0
    assert summary.namespace_counts == {}
    assert summary.classification_counts == {}
    assert summary.classification_combination_counts == {}
    assert summary.illegal_sequence_count == 0
