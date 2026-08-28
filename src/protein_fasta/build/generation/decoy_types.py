"""Runtime values and behavior required by protein-decoy generation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol

from protein_fasta.schema.build import DecoyMode

Entry = tuple[str, str]


@dataclass(frozen=True, slots=True)
class DecoyBatch:
    """Generated entries and algorithm evidence for one database build."""

    entries: tuple[Entry, ...]
    parameters: dict[str, Any]
    initial_collisions: int = 0
    unresolved_collisions: int = 0
    dropped_peptides: int = 0
    omitted_decoys: int = 0


class DecoyGeneration(Protocol):
    """Generate decoys and describe the same compiled algorithm choice."""

    @property
    def mode(self) -> DecoyMode:
        """Return the stable algorithm label."""
        ...

    @property
    def seed(self) -> int | None:
        """Return the stochastic seed, or ``None`` for deterministic reversal."""
        ...

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Generate decoys for one database's source entries."""
        ...

    def parameters(self) -> dict[str, Any]:
        """Return machine-readable provenance for this algorithm."""
        ...

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Render sentinel provenance, or the empty identity for the default."""
        ...
