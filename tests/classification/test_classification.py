"""Tests for primary and independent identifier classification."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from protein_fasta.classification.kinds import EntryKind
from protein_fasta.compile import (
    make_identifier_classifier,
    make_pattern_matcher,
    make_pattern_resolution,
)
from protein_fasta.schema.classification import (
    ExplicitPatternDocument,
    IdentifierClassificationDocument,
    IdentifierPatternDocument,
    InferredPatternDocument,
)


@pytest.mark.parametrize(
    ("identifier", "expected"),
    [
        ("aa|database", EntryKind.SENTINEL),
        ("REV_sp|P1|ONE", EntryKind.DECOY),
        ("sp|P1|ONE_p_target", EntryKind.ENTRAPMENT),
        ("sp|Cont_P1|ONE", EntryKind.CONTAMINANT),
        ("sp|P1|ONE", EntryKind.TARGET),
    ],
)
def test_primary_classifier_preserves_configured_precedence(
    identifier: str,
    expected: EntryKind,
) -> None:
    classifier = make_identifier_classifier(IdentifierClassificationDocument())
    assert classifier.classify(identifier) is expected


def test_pattern_resolution_keeps_independent_decoy_and_contaminant_matches() -> None:
    resolution = make_pattern_resolution(
        IdentifierPatternDocument(
            decoy=InferredPatternDocument(candidates=(r"^REV_", r"^DECOY_")),
            contaminant=ExplicitPatternDocument(patterns=(r"Cont_",)),
        )
    )
    resolution.observe("REV_Cont_P1")
    resolution.observe("P2")
    resolved = resolution.resolve()

    assert resolved.n_identifiers == 2
    assert resolved.decoy.patterns == (r"^REV_",)
    assert resolved.decoy.source == "inferred"
    assert [item.count for item in resolved.decoy.match_counts] == [1, 0]
    assert resolved.contaminant.patterns == (r"Cont_",)
    assert resolved.contaminant.source == "explicit"
    assert make_pattern_matcher(resolved.decoy).matches("REV_Cont_P1")
    assert make_pattern_matcher(resolved.contaminant).matches("REV_Cont_P1")
    assert resolved.model_validate_json(resolved.model_dump_json()) == resolved


def test_inferred_empty_patterns_record_none_source() -> None:
    resolved = make_pattern_resolution(
        IdentifierPatternDocument(
            decoy=InferredPatternDocument(candidates=(r"^REV_",)),
            contaminant=InferredPatternDocument(candidates=()),
        )
    ).resolve()
    assert resolved.decoy.source == "none"
    assert resolved.decoy.patterns == ()
    assert resolved.contaminant.source == "none"


def test_schema_rejects_invalid_regular_expression() -> None:
    with pytest.raises(ValidationError, match="invalid FASTA identifier regex"):
        ExplicitPatternDocument(patterns=("[",))


def test_primary_classifier_requires_a_nonempty_decoy_prefix() -> None:
    with pytest.raises(ValidationError):
        IdentifierClassificationDocument(decoy_prefix="")


def test_pattern_document_discriminator_loads_authored_data() -> None:
    document = IdentifierPatternDocument.model_validate(
        {
            "decoy": {"mode": "explicit", "patterns": ["^D_"]},
            "contaminant": {"mode": "infer", "candidates": ["^C_"]},
        }
    )
    assert isinstance(document.decoy, ExplicitPatternDocument)
    assert isinstance(document.contaminant, InferredPatternDocument)
