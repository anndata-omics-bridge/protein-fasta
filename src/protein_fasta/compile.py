"""Compile Pydantic policy documents into narrow runtime behavior."""

from __future__ import annotations

import re

from protein_fasta.classification.classifier import IdentifierClassifier
from protein_fasta.classification.patterns import (
    ExplicitPatternPolicy,
    IdentifierPatternAccumulator,
    InferredPatternPolicy,
    PatternMatcher,
    PatternPolicy,
    ResolvedPatternSet,
)
from protein_fasta.headers.generic import GenericHeaderInterpreter, HeaderInterpreter
from protein_fasta.headers.uniprot import UniProtHeaderInterpreter
from protein_fasta.schema.classification import (
    IdentifierClassificationDocument,
    IdentifierPatternDocument,
    InferredPatternDocument,
    PatternDocument,
    PatternMatchCountDocument,
    ResolvedIdentifierPatternsDocument,
    ResolvedPatternDocument,
)
from protein_fasta.schema.headers import (
    GenericHeaderDocument,
    HeaderDocument,
)
from protein_fasta.schema.validation import SequenceValidationDocument
from protein_fasta.validation.namespace import (
    IdentifierNamespaceClassifier,
    NamespaceAccumulator,
    NamespaceRule,
)
from protein_fasta.validation.sequence import ProteinSequenceValidator


def make_header_interpreter(document: HeaderDocument, /) -> HeaderInterpreter:
    """Construct the configured header interpretation."""
    if isinstance(document, GenericHeaderDocument):
        return GenericHeaderInterpreter()
    return UniProtHeaderInterpreter()


def make_identifier_classifier(
    document: IdentifierClassificationDocument,
    /,
) -> IdentifierClassifier:
    """Construct a primary identifier classifier from validated regexes."""
    return IdentifierClassifier(
        sentinel=tuple(re.compile(pattern) for pattern in document.sentinel_regexes),
        decoy_prefix=document.decoy_prefix,
        entrapment=tuple(re.compile(pattern) for pattern in document.entrapment_regexes),
        contaminant=tuple(re.compile(pattern) for pattern in document.contaminant_regexes),
    )


class PatternResolution:
    """Composition-boundary adapter producing serialized pattern provenance."""

    __slots__ = ("_accumulator",)

    def __init__(self, accumulator: IdentifierPatternAccumulator) -> None:
        """Retain the schema-free runtime accumulator."""
        self._accumulator = accumulator

    def observe(self, identifier: str, /) -> None:
        """Count candidate matches for one identifier."""
        self._accumulator.observe(identifier)

    def resolve(self) -> ResolvedIdentifierPatternsDocument:
        """Return validated, serializable resolution provenance."""
        resolved = self._accumulator.resolve()
        return ResolvedIdentifierPatternsDocument(
            decoy=self._resolved_document(resolved.decoy),
            contaminant=self._resolved_document(resolved.contaminant),
            n_identifiers=resolved.n_identifiers,
        )

    @staticmethod
    def _resolved_document(resolved: ResolvedPatternSet) -> ResolvedPatternDocument:
        return ResolvedPatternDocument(
            patterns=resolved.patterns,
            source=resolved.source.value,
            match_counts=tuple(
                PatternMatchCountDocument(pattern=item.pattern, count=item.count)
                for item in resolved.match_counts
            ),
        )


def make_pattern_resolution(document: IdentifierPatternDocument, /) -> PatternResolution:
    """Construct one streaming pattern-resolution operation."""
    accumulator = IdentifierPatternAccumulator(
        decoy=_make_pattern_policy(document.decoy),
        contaminant=_make_pattern_policy(document.contaminant),
    )
    return PatternResolution(accumulator)


def make_pattern_matcher(document: ResolvedPatternDocument, /) -> PatternMatcher:
    """Compile one resolved pattern document for repeated matching."""
    return PatternMatcher(document.patterns)


def make_sequence_validator(document: SequenceValidationDocument, /) -> ProteinSequenceValidator:
    """Construct configured protein-sequence normalization and validation."""
    return ProteinSequenceValidator(
        allowed_residues=document.standard_residues + document.tolerated_residues,
        strip_trailing_stop=document.strip_trailing_stop,
    )


def make_namespace_classifier(
    document: SequenceValidationDocument,
    /,
    *,
    decoy_prefix: str = "",
) -> IdentifierNamespaceClassifier:
    """Construct configured identifier-namespace recognition."""
    rules = tuple(
        NamespaceRule(rule.name, re.compile(rule.pattern)) for rule in document.id_namespaces
    )
    return IdentifierNamespaceClassifier(rules, decoy_prefix=decoy_prefix)


def make_namespace_accumulator(
    document: SequenceValidationDocument,
    /,
    *,
    decoy_prefix: str = "",
) -> NamespaceAccumulator:
    """Construct bounded namespace counting for one input scan."""
    return NamespaceAccumulator(
        make_namespace_classifier(document, decoy_prefix=decoy_prefix),
        max_namespaces=document.max_reported_id_namespaces,
    )


def _make_pattern_policy(document: PatternDocument) -> PatternPolicy:
    if isinstance(document, InferredPatternDocument):
        return InferredPatternPolicy(document.candidates)
    return ExplicitPatternPolicy(document.patterns)
