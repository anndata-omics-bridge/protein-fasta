"""Frozen runtime values for biological and search protein inventories."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

type BiologicalKind = Literal["sentinel", "target", "contaminant", "entrapment"]
type SourceRole = Literal["target", "contaminant", "foreign", "entrapment"]
type DecoySourceKind = Literal["target", "contaminant", "entrapment"]


@dataclass(frozen=True, slots=True)
class ProteinInventoryEntry:
    """One non-decoy biological or metadata record in final FASTA order."""

    final_order: int
    raw_header: str
    identifier: str
    description: str | None
    sequence: str
    kind: BiologicalKind
    contaminant_group: str | None
    sequence_hash: str
    entrapment_strategy: str | None
    source_order: int | None
    record_order: int | None
    source_id: str | None
    source_role: SourceRole | None


@dataclass(frozen=True, slots=True)
class DecoyInventoryEntry:
    """One generated decoy with required biological-source provenance."""

    final_order: int
    raw_header: str
    identifier: str
    description: str | None
    sequence: str
    contaminant_group: str | None
    sequence_hash: str
    entrapment_strategy: str | None
    decoy_strategy: str
    decoy_source_order: int
    decoy_source_id: str
    decoy_source_kind: DecoySourceKind

    @property
    def kind(self) -> Literal["decoy"]:
        """Return the inventory discriminator shared by search entries."""
        return "decoy"


type SearchInventoryEntry = ProteinInventoryEntry | DecoyInventoryEntry


@dataclass(frozen=True, slots=True)
class BiologicalDatabase:
    """One assembled biological database with a single canonical record tuple."""

    entries: tuple[ProteinInventoryEntry, ...]


@dataclass(frozen=True, slots=True)
class SearchDatabase:
    """One biological database plus source-linked generated decoys."""

    entries: tuple[SearchInventoryEntry, ...]
