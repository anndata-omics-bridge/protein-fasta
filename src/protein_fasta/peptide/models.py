"""Frozen runtime values for protein-to-peptide workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Protocol

from protein_fasta.analytics.digestion import Digestion

type PeptideProteinKind = Literal["target", "contaminant", "entrapment", "decoy"]


@dataclass(frozen=True, slots=True)
class PeptideProtein:
    """One classified protein supplied to theoretical digestion."""

    ordinal: int
    identifier: str
    kind: PeptideProteinKind
    sequence: str


@dataclass(frozen=True, slots=True)
class PeptideInventoryEntry:
    """One unique peptide with exact population counts."""

    peptide_id: str
    sequence: str
    length: int
    missed_cleavages: int
    mapping_count: int
    protein_count: int
    target_count: int
    contaminant_count: int
    entrapment_count: int
    decoy_count: int


@dataclass(frozen=True, slots=True)
class ProteinPeptideMapping:
    """One unique peptide-to-protein relationship."""

    peptide_id: str
    peptide_sequence: str
    protein_order: int
    protein_id: str
    protein_kind: PeptideProteinKind
    missed_cleavages: int


@dataclass(frozen=True, slots=True)
class PeptideDatabase:
    """Canonical unique peptides and their protein mappings."""

    peptides: tuple[PeptideInventoryEntry, ...]
    mappings: tuple[ProteinPeptideMapping, ...]


class PeptideExecutor(Protocol):
    """Client-owned execution capability for one theoretical digest."""

    @property
    def name(self) -> str:
        """Return the stable execution-backend name."""
        ...

    def execute(
        self,
        proteins: tuple[PeptideProtein, ...],
        digestion: Digestion,
    ) -> PeptideDatabase:
        """Build one exact peptide database."""
        ...
