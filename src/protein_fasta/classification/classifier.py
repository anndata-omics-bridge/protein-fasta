"""Configured primary classification of protein identifiers."""

from __future__ import annotations

import re

from protein_fasta.classification.kinds import EntryKind


class IdentifierClassifier:
    """Classify identifiers using one precompiled, first-match policy."""

    __slots__ = ("_contaminant", "_decoy_prefix", "_entrapment", "_sentinel")

    def __init__(
        self,
        *,
        sentinel: tuple[re.Pattern[str], ...],
        decoy_prefix: str,
        entrapment: tuple[re.Pattern[str], ...],
        contaminant: tuple[re.Pattern[str], ...],
    ) -> None:
        """Retain only compiled runtime matching values."""
        self._sentinel = sentinel
        self._decoy_prefix = decoy_prefix
        self._entrapment = entrapment
        self._contaminant = contaminant

    def classify(self, identifier: str, /) -> EntryKind:
        """Return the first matching kind under the configured precedence."""
        if any(pattern.match(identifier) for pattern in self._sentinel):
            return EntryKind.SENTINEL
        if identifier.startswith(self._decoy_prefix):
            return EntryKind.DECOY
        if any(pattern.search(identifier) for pattern in self._entrapment):
            return EntryKind.ENTRAPMENT
        if any(pattern.match(identifier) for pattern in self._contaminant):
            return EntryKind.CONTAMINANT
        return EntryKind.TARGET
