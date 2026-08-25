"""Tests for sequence and identifier validation."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from protein_fasta.compile import (
    make_namespace_accumulator,
    make_namespace_classifier,
    make_sequence_validator,
)
from protein_fasta.schema.validation import (
    IdentifierNamespaceDocument,
    SequenceValidationDocument,
)
from protein_fasta.validation.namespace import UNMATCHED_NAMESPACE
from protein_fasta.validation.sequence import NormalizedSequence


def test_sequence_normalization_and_illegal_residues() -> None:
    validator = make_sequence_validator(SequenceValidationDocument())
    normalized = validator.normalize("acdx*")
    assert normalized == NormalizedSequence("ACDX", upper_cased=True, stop_stripped=True)
    assert validator.illegal_residues(normalized.sequence) == ""
    assert validator.illegal_residues("ACD-*") == "-*"
    assert "'-'" in validator.describe_violation("-*")
    assert "translated nucleotide" in validator.describe_violation("*")


def test_sequence_validator_can_preserve_a_trailing_stop() -> None:
    validator = make_sequence_validator(SequenceValidationDocument(strip_trailing_stop=False))
    assert validator.normalize("AA*") == NormalizedSequence(
        "AA*", upper_cased=False, stop_stripped=False
    )


def test_namespace_classifier_handles_named_authority_and_unmatched_values() -> None:
    document = SequenceValidationDocument(
        id_namespaces=(IdentifierNamespaceDocument(name="custom", pattern=r"^CUST_"),)
    )
    classifier = make_namespace_classifier(document, decoy_prefix="REV_")
    assert classifier.name("REV_CUST_1") == "custom"
    assert classifier.name("pf|one") == "pf|"
    assert classifier.name("plain") == UNMATCHED_NAMESPACE


def test_namespace_accumulator_bounds_one_off_authorities() -> None:
    document = SequenceValidationDocument(id_namespaces=(), max_reported_id_namespaces=1)
    accumulator = make_namespace_accumulator(document)
    accumulator.observe("a|one")
    accumulator.observe("b|two")
    accumulator.observe("a|three")
    assert accumulator.counts() == {"a|": 2, UNMATCHED_NAMESPACE: 1}


def test_namespace_schema_rejects_bad_regex_and_nonpositive_bound() -> None:
    with pytest.raises(ValidationError, match="invalid identifier namespace regex"):
        IdentifierNamespaceDocument(name="bad", pattern="[")
    with pytest.raises(ValidationError):
        SequenceValidationDocument(max_reported_id_namespaces=0)
