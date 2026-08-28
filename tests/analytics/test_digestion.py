from __future__ import annotations

import pytest
from pydantic import ValidationError

from protein_fasta.analytics.digestion import digest_sequence, protein_residue_length
from protein_fasta.analytics_compile import make_digestion
from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.validation.sequence import normalize_sequence


def _peptides(sequence: str, config: DigestionDocument) -> list[tuple[str, int]]:
    digestion = make_digestion(config)
    return [
        (peptide.sequence, peptide.missed_cleavages)
        for peptide in digest_sequence(sequence, digestion)
    ]


def test_trypsin_cleaves_after_kr_except_before_proline() -> None:
    config = DigestionDocument(min_length=1, max_length=50, missed_cleavages=0)

    assert _peptides("AKRPQK", config) == [("AK", 0), ("RPQK", 0)]


def test_missed_cleavages_generate_zero_through_configured_maximum() -> None:
    config = DigestionDocument(min_length=1, max_length=50, missed_cleavages=2)

    assert _peptides("AKBKCK", config) == [
        ("AK", 0),
        ("AKBK", 1),
        ("AKBKCK", 2),
        ("BK", 0),
        ("BKCK", 1),
        ("CK", 0),
    ]


def test_length_bounds_are_inclusive() -> None:
    config = DigestionDocument(min_length=2, max_length=4, missed_cleavages=1)

    assert _peptides("AKBKCK", config) == [
        ("AK", 0),
        ("AKBK", 1),
        ("BK", 0),
        ("BKCK", 1),
        ("CK", 0),
    ]


def test_normalization_preserves_letters_and_splits_non_letters() -> None:
    config = DigestionDocument(min_length=1, max_length=50, missed_cleavages=0)

    normalized = normalize_sequence("ak*rpqk-Xu").sequence
    assert _peptides(normalized, config) == [
        ("AK", 0),
        ("RPQK", 0),
        ("XU", 0),
    ]
    assert protein_residue_length(normalized) == 8


def test_digestion_refuses_implicit_normalization() -> None:
    with pytest.raises(ValueError, match="explicitly normalized"):
        _peptides("ak", DigestionDocument(min_length=1))


def test_kernel_returns_repeated_candidates_for_orchestration_to_deduplicate() -> None:
    config = DigestionDocument(min_length=1, max_length=50, missed_cleavages=0)

    assert _peptides("AKAK", config) == [("AK", 0), ("AK", 0)]


@pytest.mark.parametrize(
    "values",
    [
        {"min_length": 0},
        {"min_length": 10, "max_length": 9},
        {"missed_cleavages": -1},
        {"missed_cleavages": 3},
    ],
)
def test_digestion_config_rejects_invalid_bounds(values: dict[str, int]) -> None:
    with pytest.raises(ValidationError):
        DigestionDocument.model_validate(values)
