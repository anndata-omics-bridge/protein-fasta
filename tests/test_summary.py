"""Tests for streaming FASTA summaries."""

import pytest

from protein_fasta.summary import SummaryAccumulator, summarize_sequences


def test_summary_uses_inclusive_quartiles_and_residue_counts() -> None:
    summary = summarize_sequences("A" * length for length in range(1, 9))
    assert summary.n_sequences == 8
    assert summary.length_q1 == 2.75
    assert summary.length_median == 4.5
    assert summary.length_q3 == 6.25
    assert summary.total_residues == 36
    assert summary.aa_frequencies == {"A": 36}


def test_add_length_retains_distribution_without_amino_acid_counts() -> None:
    accumulator = SummaryAccumulator()
    accumulator.add_length(3)
    summary = accumulator.summary()
    assert summary.n_sequences == 1
    assert summary.length_mean == 3
    assert summary.aa_frequencies == {}


def test_merge_and_snapshot_are_independent() -> None:
    first = SummaryAccumulator()
    first.add("A")
    snapshot = first.summary()
    second = SummaryAccumulator()
    second.add("CCC")
    first.merge(second)
    assert snapshot.n_sequences == 1
    assert first.summary().aa_frequencies == {"A": 1, "C": 3}


def test_empty_and_negative_lengths_have_explicit_results() -> None:
    accumulator = SummaryAccumulator()
    assert accumulator.summary().n_sequences == 0
    with pytest.raises(ValueError, match="nonnegative"):
        accumulator.add_length(-1)
