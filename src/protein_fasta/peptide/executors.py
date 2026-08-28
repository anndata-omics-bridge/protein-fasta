"""Selectable peptide execution behaviors with one exact result contract."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from protein_fasta.analytics.digestion import Digestion
from protein_fasta.peptide.computation import (
    DigestedProtein,
    PartitionDigestResult,
    digest_partitions,
    partition_digest_inputs,
    peptide_database_from_partitions,
)
from protein_fasta.peptide.models import (
    PeptideDatabase,
    PeptideProtein,
    PeptideProteinKind,
)
from protein_fasta.registry.backend import factory


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
class SQLitePeptideExecutor:
    """Spill partition runs through a temporary SQLite workspace."""

    workers: int
    partition_size: int
    name: str = "sqlite"

    def execute(
        self,
        proteins: tuple[PeptideProtein, ...],
        digestion: Digestion,
    ) -> PeptideDatabase:
        """Digest through the SQLite execution behavior."""
        return _execute_with_workspace(
            proteins,
            digestion,
            workers=self.workers,
            partition_size=self.partition_size,
            backend=self.name,
        )


@dataclass(frozen=True, slots=True)
class DuckDBPeptideExecutor:
    """Spill partition runs through a temporary DuckDB workspace."""

    workers: int
    partition_size: int
    name: str = "duckdb"

    def execute(
        self,
        proteins: tuple[PeptideProtein, ...],
        digestion: Digestion,
    ) -> PeptideDatabase:
        """Digest through the DuckDB execution behavior."""
        return _execute_with_workspace(
            proteins,
            digestion,
            workers=self.workers,
            partition_size=self.partition_size,
            backend=self.name,
        )


def _execute_with_workspace(
    proteins: tuple[PeptideProtein, ...],
    digestion: Digestion,
    *,
    workers: int,
    partition_size: int,
    backend: str,
) -> PeptideDatabase:
    """Persist completed partitions, then merge them in stable index order."""
    protein_by_order = {protein.ordinal: protein for protein in proteins}
    with tempfile.NamedTemporaryFile(
        prefix="protein-fasta-peptides-",
        suffix=factory.suffix_for(backend),
        delete=False,
    ) as handle:
        path = Path(handle.name)
    path.unlink(missing_ok=True)
    try:
        with factory.connect(path, backend=backend) as connection:
            connection.execute(
                "CREATE TABLE partition_rows ("
                "partition_index INTEGER NOT NULL, protein_order INTEGER NOT NULL, "
                "protein_id TEXT NOT NULL, protein_kind TEXT NOT NULL, "
                "sequence TEXT NOT NULL, missed_cleavages INTEGER NOT NULL)"
            )
            tasks = partition_digest_inputs(
                proteins,
                digestion,
                partition_size=partition_size,
            )
            for partition in digest_partitions(tasks, workers=workers):
                rows = [
                    (
                        partition.index,
                        digested.protein.ordinal,
                        digested.protein.identifier,
                        digested.protein.kind,
                        sequence,
                        missed_cleavages,
                    )
                    for digested in partition.proteins
                    for sequence, missed_cleavages in digested.peptides
                ]
                connection.executemany(
                    "INSERT INTO partition_rows VALUES (?, ?, ?, ?, ?, ?)",
                    rows,
                )
            records = connection.execute(
                "SELECT partition_index, protein_order, protein_id, protein_kind, "
                "sequence, missed_cleavages FROM partition_rows "
                "ORDER BY partition_index, protein_order, sequence"
            ).fetchall()
        partitions: dict[int, dict[tuple[int, str, str], list[tuple[str, int]]]] = {}
        for row in records:
            partition_index = int(row[0])
            protein_key = (int(row[1]), str(row[2]), str(row[3]))
            partitions.setdefault(partition_index, {}).setdefault(protein_key, []).append(
                (str(row[4]), int(row[5]))
            )
        results = tuple(
            PartitionDigestResult(
                index,
                tuple(
                    DigestedProtein(
                        PeptideProtein(
                            order,
                            protein_id,
                            cast("PeptideProteinKind", protein_kind),
                            protein_by_order[order].sequence,
                        ),
                        tuple(peptides),
                    )
                    for (order, protein_id, protein_kind), peptides in sorted(group.items())
                ),
            )
            for index, group in sorted(partitions.items())
        )
        return peptide_database_from_partitions(results)
    finally:
        path.unlink(missing_ok=True)
