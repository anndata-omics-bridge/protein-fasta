"""Compile Pydantic documents into schema-free runtime rules."""

from __future__ import annotations

import re

from protein_fasta.diagnostics.runtime import (
    DiagnosticRules,
    EntryClassifier,
    NamespaceRule,
)
from protein_fasta.schema.diagnostics import (
    DiagnosticDocument,
    EntryClassifierCatalogDocument,
)


def make_diagnostic_rules(
    document: DiagnosticDocument,
    classifiers: EntryClassifierCatalogDocument,
    /,
) -> DiagnosticRules:
    """Compile diagnostic documents for repeated scalar inspection."""
    return DiagnosticRules(
        namespace_rules=tuple(
            NamespaceRule(rule.name, re.compile(rule.pattern))
            for rule in document.identifier_namespaces
        ),
        classifiers=tuple(
            EntryClassifier(
                name=classifier.name,
                match_patterns=tuple(re.compile(pattern) for pattern in classifier.match_patterns),
                removable_prefix_patterns=tuple(
                    re.compile(pattern) for pattern in classifier.removable_prefix_patterns
                ),
                removable_suffix_patterns=tuple(
                    re.compile(pattern) for pattern in classifier.removable_suffix_patterns
                ),
            )
            for classifier in classifiers.classifiers
        ),
        allowed_residues=document.allowed_residues,
    )
