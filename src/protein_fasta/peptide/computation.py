"""Backend-free deterministic protein digestion and peptide aggregation."""

from __future__ import annotations

import hashlib
import multiprocessing
from collections import Counter
from collections.abc import Iterable, Iterator
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from dataclasses import dataclass

from protein_fasta.analytics.digestion import Digestion, digest_sequence
from protein_fasta.peptide.models import (
    PeptideDatabase,
    PeptideInventoryEntry,
    PeptideProtein,
    ProteinPeptideMapping,
)


@dataclass(frozen=True, slots=True)
class PartitionDigestTask:
    """One deterministic, independently digestible protein partition."""

    index: int
    proteins: tuple[PeptideProtein, ...]
    digestion: Digestion


@dataclass(frozen=True, slots=True)
class DigestedProtein:
    """Unique theoretical peptides observed for one protein."""

    protein: PeptideProtein
    peptides: tuple[tuple[str, int], ...]


@dataclass(frozen=True, slots=True)
class PartitionDigestResult:
    """One completed partition in source order."""

    index: int
    proteins: tuple[DigestedProtein, ...]


def peptide_id(sequence: str, /) -> str:
    """Return a stable content-derived peptide identifier."""
    digest = hashlib.blake2b(sequence.encode("ascii"), digest_size=16).hexdigest()
    return f"pep|{digest}"


def partition_digest_inputs(
    proteins: Iterable[PeptideProtein],
    digestion: Digestion,
    *,
    partition_size: int,
) -> Iterator[PartitionDigestTask]:
    """Yield coarse tasks without materializing a second protein collection."""
    if partition_size < 1:
        raise ValueError("partition_size must be at least 1")
    partition: list[PeptideProtein] = []
    index = 0
    for protein in proteins:
        partition.append(protein)
        if len(partition) != partition_size:
            continue
        yield PartitionDigestTask(index, tuple(partition), digestion)
        partition.clear()
        index += 1
    if partition:
        yield PartitionDigestTask(index, tuple(partition), digestion)


def digest_partition(task: PartitionDigestTask, /) -> PartitionDigestResult:
    """Digest one task without shared mutable state."""
    digested: list[DigestedProtein] = []
    for protein in task.proteins:
        unique: dict[str, int] = {}
        for peptide in digest_sequence(protein.sequence, task.digestion):
            previous = unique.get(peptide.sequence)
            if previous is None or peptide.missed_cleavages < previous:
                unique[peptide.sequence] = peptide.missed_cleavages
        digested.append(DigestedProtein(protein, tuple(sorted(unique.items()))))
    return PartitionDigestResult(task.index, tuple(digested))


def digest_partitions(
    tasks: Iterable[PartitionDigestTask],
    *,
    workers: int,
) -> Iterator[PartitionDigestResult]:
    """Yield results in partition order with bounded process submission."""
    if workers < 1:
        raise ValueError("workers must be at least 1")
    if workers == 1:
        yield from map(digest_partition, tasks)
        return
    task_iterator = iter(tasks)
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=multiprocessing.get_context("spawn"),
    ) as executor:
        pending: dict[Future[PartitionDigestResult], int] = {}
        ready: dict[int, PartitionDigestResult] = {}
        next_index = 0
        for _ in range(workers):
            task = next(task_iterator, None)
            if task is not None:
                pending[executor.submit(digest_partition, task)] = task.index
        while pending:
            completed, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in completed:
                del pending[future]
                result = future.result()
                ready[result.index] = result
                task = next(task_iterator, None)
                if task is not None:
                    pending[executor.submit(digest_partition, task)] = task.index
            while next_index in ready:
                yield ready.pop(next_index)
                next_index += 1


def peptide_database_from_partitions(
    partitions: Iterable[PartitionDigestResult],
) -> PeptideDatabase:
    """Merge partition results into canonical peptide and mapping tuples."""
    mappings_by_key: dict[tuple[str, str], ProteinPeptideMapping] = {}
    counts: dict[str, Counter[str]] = {}
    missed: dict[str, int] = {}
    for partition in partitions:
        for digested in partition.proteins:
            protein = digested.protein
            for sequence, missed_cleavages in digested.peptides:
                identifier = peptide_id(sequence)
                mapping = ProteinPeptideMapping(
                    peptide_id=identifier,
                    peptide_sequence=sequence,
                    protein_order=protein.ordinal,
                    protein_id=protein.identifier,
                    protein_kind=protein.kind,
                    missed_cleavages=missed_cleavages,
                )
                key = (sequence, protein.identifier)
                previous = mappings_by_key.get(key)
                if previous is None:
                    mappings_by_key[key] = mapping
                    counts.setdefault(sequence, Counter())[protein.kind] += 1
                elif (
                    previous.protein_kind != protein.kind
                    or previous.missed_cleavages != missed_cleavages
                ):
                    raise ValueError(
                        f"Mapping {key!r} has conflicting protein kind or missed-cleavage count"
                    )
                missed[sequence] = min(missed.get(sequence, missed_cleavages), missed_cleavages)
    mappings = tuple(
        sorted(
            mappings_by_key.values(),
            key=lambda row: (row.peptide_sequence, row.protein_order, row.protein_id),
        )
    )
    peptides = tuple(
        PeptideInventoryEntry(
            peptide_id=peptide_id(sequence),
            sequence=sequence,
            length=len(sequence),
            missed_cleavages=missed[sequence],
            mapping_count=sum(counts[sequence].values()),
            protein_count=sum(counts[sequence].values()),
            target_count=counts[sequence]["target"],
            contaminant_count=counts[sequence]["contaminant"],
            entrapment_count=counts[sequence]["entrapment"],
            decoy_count=counts[sequence]["decoy"],
        )
        for sequence in sorted(counts)
    )
    return PeptideDatabase(peptides, mappings)
