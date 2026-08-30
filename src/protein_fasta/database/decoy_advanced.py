"""Optional ``fdr_benchmark`` implementations of advanced decoy generation."""

from __future__ import annotations

import importlib.metadata
from dataclasses import asdict, dataclass, replace
from typing import cast

from fdr_benchmark.decoypyrat import generate_decoypyrat
from fdr_benchmark.models import (
    DecoyGenerationResult,
    DecoyPyratRequest,
    DigestionSpec,
    FastaRecord,
    ProteinShuffleRequest,
)
from fdr_benchmark.protein_shuffle import generate_shuffled_decoys

from protein_fasta.database.decoy import DecoyBatch, DecoyMode, Entry


@dataclass(frozen=True, slots=True)
class ShuffleDecoyGeneration:
    """Generate seeded whole-protein shuffles through ``fdr_benchmark``."""

    _seed: int
    _implementation_version: str

    @property
    def mode(self) -> DecoyMode:
        """Return the shuffle algorithm identity."""
        return DecoyMode.SHUFFLE

    @property
    def seed(self) -> int:
        """Return the configured shuffle seed."""
        return self._seed

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Generate one seeded shuffle for each source protein."""
        records = tuple(
            _benchmark_record(description, sequence) for description, sequence in entries
        )
        result = generate_shuffled_decoys(
            records,
            ProteinShuffleRequest(prefix=prefix, seed=self._seed),
        )
        return _batch_from_result(entries, records, result, prefix, self.parameters())

    def parameters(self) -> dict[str, object]:
        """Return shuffle-generation provenance."""
        return {
            "mode": self.mode.value,
            "seed": self._seed,
            "implementation": "fdr_benchmark",
            "implementation_version": self._implementation_version,
        }

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Render compact shuffle provenance for the FASTA sentinel."""
        del initial_collisions, dropped_peptides
        return (
            f"decoys {self.mode.value} seed {self._seed} "
            f"with fdr_benchmark {self._implementation_version}"
        )


@dataclass(frozen=True, slots=True)
class DecoyPyratGeneration:
    """Generate digestion-aware DecoyPYrat proteins through ``fdr_benchmark``."""

    _request: DecoyPyratRequest
    _implementation_version: str

    @property
    def mode(self) -> DecoyMode:
        """Return the DecoyPYrat algorithm identity."""
        return DecoyMode.DECOYPYRAT

    @property
    def seed(self) -> int:
        """Return the configured DecoyPYrat seed."""
        return self._request.seed

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Generate digestion-aware decoys and collision evidence."""
        records = tuple(
            _benchmark_record(description, sequence) for description, sequence in entries
        )
        result = generate_decoypyrat(records, replace(self._request, prefix=prefix))
        return _batch_from_result(entries, records, result, prefix, self.parameters())

    def parameters(self) -> dict[str, object]:
        """Return DecoyPYrat generation and digestion provenance."""
        return {
            "mode": self.mode.value,
            "seed": self._request.seed,
            "implementation": "fdr_benchmark",
            "implementation_version": self._implementation_version,
            "digestion": cast("dict[str, object]", asdict(self._request.digestion)),
            "max_attempts": self._request.max_attempts,
            "switch_cleavage_sites": self._request.switch_cleavage_sites,
        }

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Render compact DecoyPYrat provenance and optional collision evidence."""
        digestion = self._request.digestion
        note = (
            f"decoys {self.mode.value} seed {self._request.seed} with "
            f"fdr_benchmark {self._implementation_version} "
            f"({digestion.enzyme}, missed cleavages {digestion.missed_cleavages}, "
            f"length {digestion.min_length}-{digestion.max_length}, "
            f"max attempts {self._request.max_attempts}"
        )
        if initial_collisions is not None:
            note += f", initial collisions {initial_collisions}"
        if dropped_peptides is not None:
            note += f", dropped peptides {dropped_peptides}"
        return f"{note})"


def make_shuffle_decoy_generation(seed: int, /) -> ShuffleDecoyGeneration:
    """Construct one seeded shuffle implementation."""
    return ShuffleDecoyGeneration(seed, importlib.metadata.version("fdr_benchmark"))


def make_decoypyrat_generation(
    *,
    seed: int,
    enzyme: str,
    minimum_length: int,
    maximum_length: int,
) -> DecoyPyratGeneration:
    """Construct one digestion-aware DecoyPYrat implementation."""
    request = DecoyPyratRequest(
        digestion=DigestionSpec(
            enzyme=enzyme,
            missed_cleavages=0,
            min_length=minimum_length,
            max_length=maximum_length,
        ),
        prefix="unused",
        seed=seed,
    )
    return DecoyPyratGeneration(request, importlib.metadata.version("fdr_benchmark"))


def _batch_from_result(
    entries: tuple[Entry, ...],
    records: tuple[FastaRecord, ...],
    result: DecoyGenerationResult,
    prefix: str,
    parameters: dict[str, object],
) -> DecoyBatch:
    sequences = {pair.source_id: pair.generated_sequence for pair in result.pairs}
    generated = tuple(
        (f"{prefix}{description}", sequences[record.identifier])
        for (description, _), record in zip(entries, records, strict=True)
        if record.identifier in sequences
    )
    return DecoyBatch(
        generated,
        parameters,
        len(result.collisions.initial),
        len(result.collisions.unresolved),
        result.collisions.dropped_occurrences,
        len(result.omitted_records),
    )


def _benchmark_record(description: str, sequence: str) -> FastaRecord:
    identifier, separator, remainder = description.partition(" ")
    return FastaRecord(identifier, sequence, remainder if separator else "")
