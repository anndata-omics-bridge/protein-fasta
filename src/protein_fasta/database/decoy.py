"""Backend-independent runtime behavior for protein-decoy generation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol

type Entry = tuple[str, str]


class DecoyMode(StrEnum):
    """Stable runtime identities for supported decoy algorithms."""

    REVERSE = "reverse"
    SHUFFLE = "shuffle"
    DECOYPYRAT = "decoypyrat"


@dataclass(frozen=True, slots=True)
class DecoyBatch:
    """Generated entries and observed algorithm evidence."""

    entries: tuple[Entry, ...]
    parameters: dict[str, object]
    initial_collisions: int = 0
    unresolved_collisions: int = 0
    dropped_peptides: int = 0
    omitted_decoys: int = 0


class DecoyGeneration(Protocol):
    """Generate decoys and describe the same compiled algorithm choice."""

    @property
    def mode(self) -> DecoyMode:
        """Return the stable algorithm identity."""
        ...

    @property
    def seed(self) -> int | None:
        """Return the stochastic seed, if the algorithm uses one."""
        ...

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Generate decoys for one database's source entries."""
        ...

    def parameters(self) -> dict[str, object]:
        """Return machine-readable algorithm provenance."""
        ...

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Render sentinel provenance, or the empty default identity."""
        ...


@dataclass(frozen=True, slots=True)
class ReverseDecoyGeneration:
    """Generate deterministic whole-sequence reversals."""

    @property
    def mode(self) -> DecoyMode:
        """Return the reverse algorithm identity."""
        return DecoyMode.REVERSE

    @property
    def seed(self) -> None:
        """Return no seed because reversal is deterministic."""
        return None

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Reverse each source sequence and prefix its complete header."""
        generated = tuple(
            make_decoy(description, sequence, prefix=prefix) for description, sequence in entries
        )
        return DecoyBatch(generated, self.parameters())

    def parameters(self) -> dict[str, object]:
        """Return reverse-generation provenance."""
        return {"mode": self.mode.value}

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Return no note because reversal is the documented default."""
        del initial_collisions, dropped_peptides
        return ""


def reverse_sequence(sequence: str, /) -> str:
    """Return one whole-sequence reversal."""
    return sequence[::-1]


def make_decoy(description: str, sequence: str, *, prefix: str) -> Entry:
    """Build a reversed decoy entry from one biological entry."""
    return prefix + description, reverse_sequence(sequence)
