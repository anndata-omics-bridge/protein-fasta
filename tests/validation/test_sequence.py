"""Tests for fixed high-level sequence normalization."""

from protein_fasta.validation.sequence import NormalizedSequence, normalize_sequence


def test_normalize_sequence_uppercases_and_removes_one_terminal_stop() -> None:
    assert normalize_sequence("acdx*") == NormalizedSequence(
        "ACDX",
        upper_cased=True,
        stop_stripped=True,
    )
    assert normalize_sequence("AA**") == NormalizedSequence(
        "AA*",
        upper_cased=False,
        stop_stripped=True,
    )
    assert normalize_sequence("PEPTIDE") == NormalizedSequence(
        "PEPTIDE",
        upper_cased=False,
        stop_stripped=False,
    )
