"""Schema-free protein-FASTA diagnostic rules."""

from __future__ import annotations

import re
from dataclasses import dataclass

UNMATCHED_NAMESPACE = "unmatched"
_MAX_AUTHORITY_LENGTH = 12


@dataclass(frozen=True, slots=True)
class NamespaceRule:
    """One compiled identifier namespace."""

    name: str
    pattern: re.Pattern[str]


@dataclass(frozen=True, slots=True)
class EntryClassifier:
    """One independent label and its compiled identifier expressions."""

    name: str
    match_patterns: tuple[re.Pattern[str], ...]
    removable_prefix_patterns: tuple[re.Pattern[str], ...]
    removable_suffix_patterns: tuple[re.Pattern[str], ...]


class DiagnosticRules:
    """Inspect normalized proteins using immutable compiled rules."""

    __slots__ = ("_classifiers", "_delete_allowed", "_namespace_rules")

    def __init__(
        self,
        *,
        namespace_rules: tuple[NamespaceRule, ...],
        classifiers: tuple[EntryClassifier, ...],
        allowed_residues: str,
    ) -> None:
        """Retain compiled identifier rules and the allowed-residue table."""
        self._namespace_rules = namespace_rules
        self._classifiers = classifiers
        self._delete_allowed = str.maketrans("", "", allowed_residues)

    def diagnose_identifier(self, identifier: str, /) -> tuple[str, frozenset[str]]:
        """Return the undecorated namespace and every matching classification."""
        working, classifications = self._peel(identifier)
        for classifier in self._classifiers:
            if any(pattern.search(working) for pattern in classifier.match_patterns):
                classifications.add(classifier.name)
        return self._namespace(working), frozenset(classifications)

    def illegal_residues(self, sequence: str, /) -> str:
        """Return every residue outside the configured alphabet."""
        return sequence.translate(self._delete_allowed)

    def _peel(self, identifier: str) -> tuple[str, set[str]]:
        working = identifier
        classifications: set[str] = set()
        while working:
            before_pass = working
            for classifier in self._classifiers:
                for pattern in classifier.removable_prefix_patterns:
                    match = pattern.match(working)
                    if match is not None:
                        working = working[match.end() :]
                        classifications.add(classifier.name)
                for pattern in classifier.removable_suffix_patterns:
                    match = pattern.search(working)
                    if match is not None:
                        working = working[: match.start()]
                        classifications.add(classifier.name)
            if working == before_pass:
                break
        return working, classifications

    def _namespace(self, identifier: str) -> str:
        for rule in self._namespace_rules:
            if rule.pattern.match(identifier):
                return rule.name
        authority, separator, _ = identifier.partition("|")
        if (
            separator
            and authority
            and len(authority) <= _MAX_AUTHORITY_LENGTH
            and authority.isalnum()
        ):
            return f"{authority}|"
        return UNMATCHED_NAMESPACE
