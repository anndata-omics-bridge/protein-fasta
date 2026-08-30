"""Dependency-free runtime types for biological entrapment generation."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol

type Entry = tuple[str, str]


class EntrapmentStrategy(StrEnum):
    """Stable runtime identities for biological entrapment algorithms."""

    SHUFFLED = "shuffled"
    FOREIGN_SPECIES = "foreign-species"


@dataclass(frozen=True, slots=True)
class EntrapmentPeptidePair:
    """One target-to-entrapment peptide relationship retained as evidence."""

    source_id: str
    target_peptide: str
    generated_peptide: str
    fold_index: int


@dataclass(frozen=True, slots=True)
class EntrapmentBatch:
    """Generated entries, peptide relationships, and algorithm evidence."""

    entries: tuple[Entry, ...]
    peptide_pairs: tuple[EntrapmentPeptidePair, ...]
    parameters: dict[str, object]
    requested_fold: int
    achieved_fold: int
    failures: int
    proteins_affected: int
    source_proteins: int

    @property
    def is_complete(self) -> bool:
        """Return whether every source protein reached the requested fold."""
        return not self.failures and self.achieved_fold >= self.requested_fold

    @property
    def complete_proteins(self) -> int:
        """Return the number of proteins that reached the requested fold."""
        return self.source_proteins - self.proteins_affected


class EntrapmentGeneration(Protocol):
    """Generate and describe one compiled biological-entrapment algorithm."""

    @property
    def strategy(self) -> EntrapmentStrategy:
        """Return the stable strategy identity."""
        ...

    @property
    def seed(self) -> int:
        """Return the stochastic seed."""
        ...

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]:
        """Apply the sequence policy required by the generated target space."""
        ...

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch:
        """Generate entrapment entries from one biological source space."""
        ...

    def parameters(self) -> dict[str, object]:
        """Return machine-readable strategy provenance."""
        ...

    def annotation(self, batch: EntrapmentBatch | None = None) -> str:
        """Render sentinel provenance for this strategy and observed batch."""
        ...


def format_entrapment_peptide_pairs(pairs: Iterable[EntrapmentPeptidePair], /) -> str:
    """Serialize entrapment peptide relationships as stable TSV."""
    lines = ["source_id\tfold_index\ttarget_peptide\tgenerated_peptide"]
    lines.extend(
        "\t".join(
            (
                pair.source_id,
                str(pair.fold_index),
                pair.target_peptide,
                pair.generated_peptide,
            )
        )
        for pair in pairs
    )
    return "\n".join(lines) + "\n"
