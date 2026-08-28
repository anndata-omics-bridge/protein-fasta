"""Tests for compiled scalar FASTA diagnostics."""

from protein_fasta.compile import make_diagnostic_rules
from protein_fasta.documents import (
    load_builtin_diagnostic_document,
    load_builtin_entry_classifier_document,
)


def test_diagnostics_peel_stacked_decorations_and_keep_overlapping_labels() -> None:
    rules = make_diagnostic_rules(
        load_builtin_diagnostic_document(),
        load_builtin_entry_classifier_document(),
    )

    namespace, labels = rules.diagnose_identifier("REV_CON__sp|P12345|RL40_YEAST")

    assert namespace == "uniprot"
    assert labels == frozenset({"decoy", "contaminant"})


def test_diagnostics_match_undecorated_identifier_and_report_illegal_residues() -> None:
    rules = make_diagnostic_rules(
        load_builtin_diagnostic_document(),
        load_builtin_entry_classifier_document(),
    )

    namespace, labels = rules.diagnose_identifier("REV_CON__sp|P1|ENTRY_ORG")

    assert namespace == "uniprot"
    assert labels == frozenset({"decoy", "contaminant"})
    assert rules.illegal_residues("ACD-*?") == "-*?"
