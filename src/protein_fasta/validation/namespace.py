"""Configured recognition and bounded counting of identifier namespaces."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass

UNMATCHED_NAMESPACE = "unmatched"
_MAX_AUTHORITY_LENGTH = 12


@dataclass(frozen=True, slots=True)
class NamespaceRule:
    """One compiled, named identifier namespace."""

    name: str
    pattern: re.Pattern[str]


class IdentifierNamespaceClassifier:
    """Name the configured namespace of protein identifier tokens."""

    __slots__ = ("_decoy_prefix", "_rules")

    def __init__(self, rules: tuple[NamespaceRule, ...], *, decoy_prefix: str = "") -> None:
        """Retain rules in precedence order."""
        self._rules = rules
        self._decoy_prefix = decoy_prefix

    def name(self, identifier: str, /) -> str:
        """Return a named namespace, authority prefix, or `unmatched`."""
        token = identifier.removeprefix(self._decoy_prefix) if self._decoy_prefix else identifier
        for rule in self._rules:
            if rule.pattern.match(token):
                return rule.name
        authority, separator, _ = token.partition("|")
        if (
            separator
            and authority
            and len(authority) <= _MAX_AUTHORITY_LENGTH
            and authority.isalnum()
        ):
            return f"{authority}|"
        return UNMATCHED_NAMESPACE


class NamespaceAccumulator:
    """Count namespaces while bounding one-off authority rows."""

    __slots__ = ("_classifier", "_counts", "_max_namespaces")

    def __init__(
        self,
        classifier: IdentifierNamespaceClassifier,
        *,
        max_namespaces: int,
    ) -> None:
        """Configure namespace interpretation and output cardinality."""
        self._classifier = classifier
        self._max_namespaces = max_namespaces
        self._counts: Counter[str] = Counter()

    def observe(self, identifier: str, /) -> None:
        """Count one identifier without exceeding the configured row bound."""
        namespace = self._classifier.name(identifier)
        if namespace not in self._counts and len(self._counts) >= self._max_namespaces:
            self._counts[UNMATCHED_NAMESPACE] += 1
            return
        self._counts[namespace] += 1

    def counts(self) -> dict[str, int]:
        """Return an independent count snapshot."""
        return dict(self._counts)
