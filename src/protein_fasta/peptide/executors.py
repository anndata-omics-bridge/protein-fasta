"""Selectable peptide execution behaviors with one exact result contract."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Protocol

from protein_fasta.analytics.digestion import Digestion
from protein_fasta.peptide.computation import (
    PartitionDigestResult,
    digest_partitions,
    partition_digest_inputs,
    peptide_database_from_partitions,
)
from protein_fasta.peptide.models import PeptideDatabase, PeptideProtein


class PartitionWorkspace(Protocol):
    """Persist and merge already-computed peptide partitions."""

    def merge(
        self,
        partitions: Iterable[PartitionDigestResult],
        proteins: tuple[PeptideProtein, ...],
    ) -> PeptideDatabase:
        """Return one peptide database from partitions in stable index order."""
        ...


@dataclass(frozen=True, slots=True)
class MemoryPeptideExecutor:
    """Merge partition results directly in memory."""

    workers: int
    partition_size: int
    name: str = "memory"

    def execute(
        self,
        proteins: tuple[PeptideProtein, ...],
        digestion: Digestion,
    ) -> PeptideDatabase:
        """Digest and merge deterministic in-memory partitions."""
        tasks = partition_digest_inputs(
            proteins,
            digestion,
            partition_size=self.partition_size,
        )
        return peptide_database_from_partitions(digest_partitions(tasks, workers=self.workers))


@dataclass(frozen=True, slots=True)
class WorkspacePeptideExecutor:
    """Digest partitions and merge them through an injected temporary workspace."""

    workers: int
    partition_size: int
    name: str
    workspace: PartitionWorkspace

    def execute(
        self,
        proteins: tuple[PeptideProtein, ...],
        digestion: Digestion,
    ) -> PeptideDatabase:
        """Digest and hand completed partitions to the workspace."""
        tasks = partition_digest_inputs(
            proteins,
            digestion,
            partition_size=self.partition_size,
        )
        return self.workspace.merge(
            digest_partitions(tasks, workers=self.workers),
            proteins,
        )
