"""Native shuffled and DecoyPYrat protein-decoy generation.

Derived from FDR Benchmark revision ``bbf582e`` and modified for the native
``protein_fasta`` runtime model and digestion implementation.
"""

from __future__ import annotations

import importlib.metadata
import random
from dataclasses import dataclass

from protein_fasta.analytics.digestion import (
    Digestion,
    digest_segments,
    digest_sequence,
    is_selected_peptide,
    peptide_universe,
)
from protein_fasta.database.collisions import shuffled_candidate
from protein_fasta.database.decoy import DecoyBatch, DecoyMode, Entry

_MAX_ATTEMPTS = 100
_CLEAVAGE_RESIDUES = ("K", "R")
_PRECEDING_PROBE = "AK"
_FOLLOWING_PROBE = "PA"
_IMPLEMENTATION_VERSION = importlib.metadata.version("protein-fasta")


@dataclass(frozen=True, slots=True)
class ShuffleDecoyGeneration:
    """Generate seeded whole-protein shuffled decoys."""

    _seed: int

    @property
    def mode(self) -> DecoyMode:
        """Return the shuffle algorithm identity."""
        return DecoyMode.SHUFFLE

    @property
    def seed(self) -> int:
        """Return the configured seed."""
        return self._seed

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Shuffle every complete protein with an isolated random generator."""
        _require_unique_ids(entries)
        rng = random.Random(self._seed)
        generated: list[Entry] = []
        for description, sequence in entries:
            residues = list(sequence)
            rng.shuffle(residues)
            generated.append((prefix + description, "".join(residues)))
        return DecoyBatch(tuple(generated), self.parameters())

    def parameters(self) -> dict[str, object]:
        """Return shuffled-decoy provenance."""
        return {
            "mode": self.mode.value,
            "seed": self._seed,
            "implementation": "protein_fasta",
            "implementation_version": _IMPLEMENTATION_VERSION,
        }

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Render compact shuffled-decoy provenance."""
        del initial_collisions, dropped_peptides
        return (
            f"decoys {self.mode.value} seed {self._seed} "
            f"with protein_fasta {_IMPLEMENTATION_VERSION}"
        )


@dataclass(frozen=True, slots=True)
class DecoyPyratGeneration:
    """Generate digestion-aware DecoyPYrat proteins."""

    _seed: int
    _digestion: Digestion

    @property
    def mode(self) -> DecoyMode:
        """Return the DecoyPYrat algorithm identity."""
        return DecoyMode.DECOYPYRAT

    @property
    def seed(self) -> int:
        """Return the configured seed."""
        return self._seed

    def generate(self, entries: tuple[Entry, ...], *, prefix: str) -> DecoyBatch:
        """Generate DecoyPYrat decoys and collision evidence."""
        _require_unique_ids(entries)
        target_peptides = peptide_universe(
            (sequence for _, sequence in entries),
            self._digestion,
        )
        initial_sequences = tuple(reverse_and_switch(sequence) for _, sequence in entries)
        initial_collisions = {
            peptide.sequence
            for sequence in initial_sequences
            for peptide in digest_sequence(sequence, self._digestion)
            if peptide.sequence in target_peptides
        }
        resolver = _CollisionResolver(
            target_peptides,
            random.Random(self._seed),
            self._digestion,
        )
        dropped_occurrences = 0
        omitted_decoys = 0
        generated: list[Entry] = []
        for (description, _), initial_sequence in zip(
            entries,
            initial_sequences,
            strict=True,
        ):
            corrected_segments: list[str] = []
            for peptide in digest_segments(initial_sequence, self._digestion):
                if not (
                    is_selected_peptide(peptide, self._digestion) and peptide in target_peptides
                ):
                    corrected_segments.append(peptide)
                    continue
                corrected = resolver.resolve(peptide)
                if corrected is None:
                    dropped_occurrences += 1
                else:
                    corrected_segments.append(corrected)
            decoy_sequence = "".join(corrected_segments)
            if not decoy_sequence:
                omitted_decoys += 1
                continue
            generated.append((prefix + description, decoy_sequence))

        final_collisions = {
            peptide.sequence
            for _, sequence in generated
            for peptide in digest_sequence(sequence, self._digestion)
            if peptide.sequence in target_peptides
        }
        return DecoyBatch(
            tuple(generated),
            self.parameters(),
            initial_collisions=len(initial_collisions),
            unresolved_collisions=len(final_collisions),
            dropped_peptides=dropped_occurrences,
            omitted_decoys=omitted_decoys,
        )

    def parameters(self) -> dict[str, object]:
        """Return DecoyPYrat generation and digestion provenance."""
        return {
            "mode": self.mode.value,
            "seed": self._seed,
            "implementation": "protein_fasta",
            "implementation_version": _IMPLEMENTATION_VERSION,
            "digestion": _digestion_parameters(self._digestion),
            "max_attempts": _MAX_ATTEMPTS,
            "switch_cleavage_sites": True,
        }

    def annotation(
        self,
        *,
        initial_collisions: int | None = None,
        dropped_peptides: int | None = None,
    ) -> str:
        """Render DecoyPYrat provenance and observed collision evidence."""
        note = (
            f"decoys {self.mode.value} seed {self._seed} with "
            f"protein_fasta {_IMPLEMENTATION_VERSION} "
            f"({self._digestion.cleavage.pattern}, missed cleavages 0, "
            f"length {self._digestion.min_length}-{self._digestion.max_length}, "
            f"max attempts {_MAX_ATTEMPTS}"
        )
        if initial_collisions is not None:
            note += f", initial collisions {initial_collisions}"
        if dropped_peptides is not None:
            note += f", dropped peptides {dropped_peptides}"
        return f"{note})"


def make_shuffle_decoy_generation(seed: int, /) -> ShuffleDecoyGeneration:
    """Construct one seeded whole-protein shuffle implementation."""
    return ShuffleDecoyGeneration(seed)


def make_decoypyrat_generation(
    *,
    seed: int,
    digestion: Digestion,
) -> DecoyPyratGeneration:
    """Construct one digestion-aware DecoyPYrat implementation."""
    return DecoyPyratGeneration(
        seed,
        Digestion(
            enzyme=digestion.enzyme,
            cleavage=digestion.cleavage,
            min_length=digestion.min_length,
            max_length=digestion.max_length,
            missed_cleavages=0,
        ),
    )


def reverse_and_switch(
    sequence: str,
    *,
    cleavage_residues: tuple[str, ...] = _CLEAVAGE_RESIDUES,
    switch: bool = True,
) -> str:
    """Reverse a protein and move reversed cleavage residues one place left."""
    reversed_sequence = list(sequence[::-1])
    if not switch:
        return "".join(reversed_sequence)
    sites = frozenset(cleavage_residues)
    for index, residue in enumerate(reversed_sequence):
        if index and residue in sites:
            reversed_sequence[index - 1], reversed_sequence[index] = (
                residue,
                reversed_sequence[index - 1],
            )
    return "".join(reversed_sequence)


class _CollisionResolver:
    """Resolve each distinct target/decoy peptide collision once."""

    def __init__(
        self,
        target_peptides: frozenset[str],
        rng: random.Random,
        digestion: Digestion,
    ) -> None:
        self._targets = target_peptides
        self._rng = rng
        self._digestion = digestion
        self._replacements: dict[str, str] = {}
        self._impossible: set[str] = set()

    def resolve(self, peptide: str) -> str | None:
        """Return an admissible replacement, or no value when none exists."""
        if peptide in self._replacements:
            return self._replacements[peptide]
        if peptide in self._impossible:
            return None
        expected = self._cleavage_signature(peptide)
        candidate, _ = shuffled_candidate(
            peptide,
            rng=self._rng,
            forbidden=self._targets,
            fix_c_term=True,
            accepts=lambda value: self._cleavage_signature(value) == expected,
            max_attempts=_MAX_ATTEMPTS,
        )
        if candidate is None:
            self._impossible.add(peptide)
            return None
        self._replacements[peptide] = candidate
        return candidate

    def _cleavage_signature(self, value: str) -> tuple[int, int, int]:
        """Describe cleavage alone and against preceding/following probes."""
        return (
            len(digest_segments(value, self._digestion)),
            len(digest_segments(_PRECEDING_PROBE + value, self._digestion)),
            len(digest_segments(value + _FOLLOWING_PROBE, self._digestion)),
        )


def _digestion_parameters(digestion: Digestion) -> dict[str, object]:
    return {
        "enzyme": digestion.cleavage.pattern,
        "missed_cleavages": digestion.missed_cleavages,
        "min_length": digestion.min_length,
        "max_length": digestion.max_length,
        "normalize_i_to_l": False,
        "cleavage_residues": list(_CLEAVAGE_RESIDUES),
    }


def _require_unique_ids(entries: tuple[Entry, ...]) -> None:
    identifiers = [_identifier(description) for description, _ in entries]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Input FASTA identifiers must be unique.")


def _identifier(description: str) -> str:
    return description.partition(" ")[0]
