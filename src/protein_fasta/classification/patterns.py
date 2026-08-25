"""Resolve inferred or explicit identifier-pattern policies."""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol


class PatternSource(StrEnum):
    """How one effective identifier pattern set was obtained."""

    EXPLICIT = "explicit"
    INFERRED = "inferred"
    NONE = "none"


@dataclass(frozen=True, slots=True)
class PatternMatchCount:
    """Observed match count for one tested regular expression."""

    pattern: str
    count: int


@dataclass(frozen=True, slots=True)
class ResolvedPatternSet:
    """Effective patterns and the complete candidate-count audit."""

    patterns: tuple[str, ...]
    source: PatternSource
    match_counts: tuple[PatternMatchCount, ...]


@dataclass(frozen=True, slots=True)
class ResolvedIdentifierPatterns:
    """Resolved decoy and contaminant patterns for one FASTA input."""

    decoy: ResolvedPatternSet
    contaminant: ResolvedPatternSet
    n_identifiers: int


class PatternPolicy(Protocol):
    """Resolve observed candidate counts to one effective pattern set."""

    @property
    def tested_patterns(self) -> tuple[str, ...]:
        """Return every pattern whose matches must be counted."""
        ...

    def resolve(self, counts: Mapping[str, int], /) -> ResolvedPatternSet:
        """Return an effective immutable pattern set."""
        ...


@dataclass(frozen=True, slots=True)
class InferredPatternPolicy:
    """Retain exactly those candidate patterns observed at least once."""

    candidates: tuple[str, ...]

    @property
    def tested_patterns(self) -> tuple[str, ...]:
        """Return candidate patterns in authored order."""
        return self.candidates

    def resolve(self, counts: Mapping[str, int], /) -> ResolvedPatternSet:
        """Retain candidates that matched at least one identifier."""
        patterns = tuple(pattern for pattern in self.candidates if counts[pattern])
        return ResolvedPatternSet(
            patterns=patterns,
            source=PatternSource.INFERRED if patterns else PatternSource.NONE,
            match_counts=_ordered_counts(self.candidates, counts),
        )


@dataclass(frozen=True, slots=True)
class ExplicitPatternPolicy:
    """Retain exactly the authored patterns, including an empty set."""

    patterns: tuple[str, ...]

    @property
    def tested_patterns(self) -> tuple[str, ...]:
        """Return explicit patterns in authored order."""
        return self.patterns

    def resolve(self, counts: Mapping[str, int], /) -> ResolvedPatternSet:
        """Return the explicit patterns and their observed counts."""
        return ResolvedPatternSet(
            patterns=self.patterns,
            source=PatternSource.EXPLICIT,
            match_counts=_ordered_counts(self.patterns, counts),
        )


class IdentifierPatternAccumulator:
    """Count candidate matches while raw protein identifiers stream past."""

    __slots__ = (
        "_contaminant",
        "_contaminant_counts",
        "_contaminant_regexes",
        "_decoy",
        "_decoy_counts",
        "_decoy_regexes",
        "_n_identifiers",
    )

    def __init__(self, decoy: PatternPolicy, contaminant: PatternPolicy) -> None:
        """Compile each tested regular expression exactly once."""
        self._decoy = decoy
        self._contaminant = contaminant
        self._decoy_counts = dict.fromkeys(decoy.tested_patterns, 0)
        self._contaminant_counts = dict.fromkeys(contaminant.tested_patterns, 0)
        self._decoy_regexes = tuple(re.compile(pattern) for pattern in decoy.tested_patterns)
        self._contaminant_regexes = tuple(
            re.compile(pattern) for pattern in contaminant.tested_patterns
        )
        self._n_identifiers = 0

    def observe(self, identifier: str, /) -> None:
        """Count every configured pattern matching one raw identifier."""
        self._n_identifiers += 1
        self._observe(
            identifier, self._decoy.tested_patterns, self._decoy_regexes, self._decoy_counts
        )
        self._observe(
            identifier,
            self._contaminant.tested_patterns,
            self._contaminant_regexes,
            self._contaminant_counts,
        )

    def resolve(self) -> ResolvedIdentifierPatterns:
        """Return the effective patterns and complete count audit."""
        return ResolvedIdentifierPatterns(
            decoy=self._decoy.resolve(self._decoy_counts),
            contaminant=self._contaminant.resolve(self._contaminant_counts),
            n_identifiers=self._n_identifiers,
        )

    @staticmethod
    def _observe(
        identifier: str,
        patterns: tuple[str, ...],
        regexes: tuple[re.Pattern[str], ...],
        counts: dict[str, int],
    ) -> None:
        for pattern, regex in zip(patterns, regexes, strict=True):
            if regex.search(identifier):
                counts[pattern] += 1


class PatternMatcher:
    """Match values against one immutable compiled pattern set."""

    __slots__ = ("_patterns",)

    def __init__(self, patterns: tuple[str, ...]) -> None:
        """Compile patterns once for repeated matching."""
        self._patterns = tuple(re.compile(pattern) for pattern in patterns)

    def matches(self, value: str, /) -> bool:
        """Return whether at least one configured pattern matches."""
        return any(pattern.search(value) for pattern in self._patterns)


def _ordered_counts(
    patterns: tuple[str, ...],
    counts: Mapping[str, int],
) -> tuple[PatternMatchCount, ...]:
    return tuple(PatternMatchCount(pattern, counts[pattern]) for pattern in patterns)
