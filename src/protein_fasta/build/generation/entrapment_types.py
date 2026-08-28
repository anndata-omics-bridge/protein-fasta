"""Dependency-free runtime values exchanged with entrapment implementations."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Protocol

from protein_fasta.schema.build import EntrapmentStrategy

Entry = tuple[str, str]


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
    parameters: dict[str, Any]
    requested_fold: int
    achieved_fold: int
    failures: int
    proteins_affected: int
    source_proteins: int

    @property
    def is_complete(self) -> bool:
        """Whether every source protein reached the requested fold."""
        return not self.failures and self.achieved_fold >= self.requested_fold

    @property
    def complete_proteins(self) -> int:
        """Source proteins that reached the requested fold in full."""
        return self.source_proteins - self.proteins_affected


class EntrapmentGeneration(Protocol):
    """Generate and describe one compiled entrapment algorithm."""

    @property
    def strategy(self) -> EntrapmentStrategy: ...

    @property
    def seed(self) -> int: ...

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]: ...

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch: ...

    def parameters(self) -> dict[str, Any]: ...

    def annotation(self, batch: EntrapmentBatch | None = None) -> str: ...


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
