"""Native shuffled-peptide and foreign-species entrapment generation.

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
from protein_fasta.database.entrapment import (
    EntrapmentBatch,
    EntrapmentPeptidePair,
    EntrapmentStrategy,
    Entry,
)

_MAX_ATTEMPTS = 100
_ENTRAPMENT_SUFFIX = "_p_target"
_IMPLEMENTATION_VERSION = importlib.metadata.version("protein-fasta")


@dataclass(frozen=True, slots=True)
class _PreparedProtein:
    description: str
    identifier: str
    segments: tuple[str, ...]
    candidates: tuple[tuple[str, ...], ...]
    achieved_fold: int


@dataclass(frozen=True, slots=True)
class ShuffledEntrapmentGeneration:
    """Generate paired entrapment proteins from constrained peptide shuffles."""

    _fold: int
    _seed: int
    _digestion: Digestion
    _normalize_i_to_l: bool
    _fix_peptide_n_term: bool
    _fix_peptide_c_term: bool

    @property
    def strategy(self) -> EntrapmentStrategy:
        """Return the shuffled strategy identity."""
        return EntrapmentStrategy.SHUFFLED

    @property
    def seed(self) -> int:
        """Return the configured seed."""
        return self._seed

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]:
        """Apply the configured I/L normalization to biological entries."""
        if not self._normalize_i_to_l:
            return entries
        return tuple((description, sequence.replace("I", "L")) for description, sequence in entries)

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch:
        """Generate paired entrapment proteins; foreign entries are unused."""
        del foreign_entries
        _require_unique_ids(entries)
        target_peptides = peptide_universe(
            (sequence for _, sequence in entries),
            self._digestion,
        )
        rng = random.Random(self._seed)
        cache: dict[str, tuple[str, ...]] = {}
        reserved: set[str] = set()
        failures: list[str] = []
        prepared: list[_PreparedProtein] = []
        for description, sequence in entries:
            protein, protein_failures = self._prepare_protein(
                description,
                sequence,
                target_peptides,
                cache,
                reserved,
                rng,
            )
            failures.extend(protein_failures)
            prepared.append(protein)
        output, peptide_pairs = self._reconstruct(prepared)
        achieved = min((protein.achieved_fold for protein in prepared), default=0)
        return EntrapmentBatch(
            tuple(output),
            tuple(peptide_pairs),
            self.parameters(),
            self._fold,
            achieved,
            len(failures),
            len(set(failures)),
            len(entries),
        )

    def parameters(self) -> dict[str, object]:
        """Return shuffled-entrapment provenance."""
        return {
            "strategy": self.strategy.value,
            "fold": self._fold,
            "seed": self._seed,
            "suffix": _ENTRAPMENT_SUFFIX,
            "digestion": _digestion_parameters(
                self._digestion,
                normalize_i_to_l=self._normalize_i_to_l,
            ),
            "fix_peptide_n_term": self._fix_peptide_n_term,
            "fix_peptide_c_term": self._fix_peptide_c_term,
            "reject_shared_foreign_proteins": False,
            "exhaustion_policy": "emit-partial",
            "implementation": "protein_fasta",
            "implementation_version": _IMPLEMENTATION_VERSION,
        }

    def annotation(self, batch: EntrapmentBatch | None = None) -> str:
        """Render shuffled-entrapment provenance and observed completeness."""
        return _annotation(
            self.strategy,
            self._fold,
            self._seed,
            self._digestion,
            self._normalize_i_to_l,
            batch,
        )

    def _prepare_protein(
        self,
        description: str,
        sequence: str,
        target_peptides: frozenset[str],
        cache: dict[str, tuple[str, ...]],
        reserved: set[str],
        rng: random.Random,
    ) -> tuple[_PreparedProtein, list[str]]:
        identifier = _identifier(description)
        segments = digest_segments(sequence, self._digestion)
        candidate_groups: list[tuple[str, ...]] = []
        failures: list[str] = []
        achieved_fold = self._fold
        for peptide in segments:
            if not is_selected_peptide(peptide, self._digestion):
                candidate_groups.append((peptide,) * self._fold)
                continue
            if peptide not in cache:
                cache[peptide] = self._generate_candidates(
                    peptide,
                    target_peptides,
                    reserved,
                    rng,
                )
            candidates = cache[peptide]
            candidate_groups.append(candidates)
            achieved_fold = min(achieved_fold, len(candidates))
            if len(candidates) < self._fold:
                failures.append(identifier)
        return (
            _PreparedProtein(
                description,
                identifier,
                segments,
                tuple(candidate_groups),
                achieved_fold,
            ),
            failures,
        )

    def _generate_candidates(
        self,
        peptide: str,
        target_peptides: frozenset[str],
        reserved: set[str],
        rng: random.Random,
    ) -> tuple[str, ...]:
        generated: list[str] = []
        for _ in range(self._fold):
            candidate, _ = shuffled_candidate(
                peptide,
                rng=rng,
                forbidden=target_peptides,
                reserved=reserved,
                fix_n_term=self._fix_peptide_n_term,
                fix_c_term=self._fix_peptide_c_term,
                max_attempts=_MAX_ATTEMPTS,
            )
            if candidate is None:
                break
            generated.append(candidate)
            reserved.add(candidate)
        return tuple(generated)

    def _reconstruct(
        self,
        prepared: list[_PreparedProtein],
    ) -> tuple[list[Entry], list[EntrapmentPeptidePair]]:
        output: list[Entry] = []
        peptide_pairs: list[EntrapmentPeptidePair] = []
        for protein in prepared:
            for fold_index in range(protein.achieved_fold):
                generated_sequence = "".join(
                    candidates[fold_index] for candidates in protein.candidates
                )
                generated_id = _entrapment_id(
                    protein.identifier,
                    fold_index,
                    self._fold,
                )
                output.append(
                    (
                        f"{generated_id} entrapment of {protein.description}",
                        generated_sequence,
                    )
                )
                peptide_pairs.extend(
                    EntrapmentPeptidePair(
                        source_id=protein.identifier,
                        target_peptide=peptide,
                        generated_peptide=candidates[fold_index],
                        fold_index=fold_index,
                    )
                    for peptide, candidates in zip(
                        protein.segments,
                        protein.candidates,
                        strict=True,
                    )
                )
        return output, peptide_pairs


@dataclass(frozen=True, slots=True)
class ForeignSpeciesEntrapmentGeneration:
    """Select entrapment proteins from a supplied foreign-species database."""

    _fold: int
    _seed: int
    _digestion: Digestion
    _normalize_i_to_l: bool
    _reject_shared_foreign: bool

    @property
    def strategy(self) -> EntrapmentStrategy:
        """Return the foreign-species strategy identity."""
        return EntrapmentStrategy.FOREIGN_SPECIES

    @property
    def seed(self) -> int:
        """Return the configured seed."""
        return self._seed

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]:
        """Apply the configured I/L normalization to biological entries."""
        if not self._normalize_i_to_l:
            return entries
        return tuple((description, sequence.replace("I", "L")) for description, sequence in entries)

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch:
        """Select deterministic foreign entrapment proteins."""
        if not foreign_entries:
            raise ValueError(
                "Foreign-species entrapment needs a source database to draw proteins from"
            )
        _require_unique_ids(entries + foreign_entries)
        target_peptides = peptide_universe(
            (sequence for _, sequence in entries),
            self._digestion,
        )
        eligible: list[Entry] = []
        failures: list[str] = []
        for description, sequence in foreign_entries:
            shared = target_peptides.intersection(
                peptide.sequence for peptide in digest_sequence(sequence, self._digestion)
            )
            if shared and self._reject_shared_foreign:
                failures.append(_identifier(description))
            else:
                eligible.append((description, sequence))
        rng = random.Random(self._seed)
        rng.shuffle(eligible)
        selected = eligible[: self._fold * len(entries)]
        output = tuple(
            (
                f"{_identifier(description)}{_ENTRAPMENT_SUFFIX} "
                f"foreign-species entrapment; source={description}",
                sequence,
            )
            for description, sequence in selected
        )
        achieved = len(output) // len(entries) if entries else 0
        return EntrapmentBatch(
            output,
            (),
            self.parameters(),
            self._fold,
            achieved,
            len(failures),
            len(set(failures)),
            len(entries),
        )

    def parameters(self) -> dict[str, object]:
        """Return foreign-species entrapment provenance."""
        return {
            "strategy": self.strategy.value,
            "fold": self._fold,
            "seed": self._seed,
            "suffix": _ENTRAPMENT_SUFFIX,
            "digestion": _digestion_parameters(
                self._digestion,
                normalize_i_to_l=self._normalize_i_to_l,
            ),
            "fix_peptide_n_term": True,
            "fix_peptide_c_term": True,
            "reject_shared_foreign_proteins": self._reject_shared_foreign,
            "exhaustion_policy": "emit-partial",
            "implementation": "protein_fasta",
            "implementation_version": _IMPLEMENTATION_VERSION,
        }

    def annotation(self, batch: EntrapmentBatch | None = None) -> str:
        """Render foreign-species provenance and observed completeness."""
        return _annotation(
            self.strategy,
            self._fold,
            self._seed,
            self._digestion,
            self._normalize_i_to_l,
            batch,
        )


def make_shuffled_entrapment_generation(
    *,
    fold: int,
    seed: int,
    digestion: Digestion,
    normalize_i_to_l: bool,
    fix_peptide_n_term: bool,
    fix_peptide_c_term: bool,
) -> ShuffledEntrapmentGeneration:
    """Construct peptide-shuffled entrapment behavior."""
    return ShuffledEntrapmentGeneration(
        fold,
        seed,
        digestion,
        normalize_i_to_l,
        fix_peptide_n_term,
        fix_peptide_c_term,
    )


def make_foreign_species_entrapment_generation(
    *,
    fold: int,
    seed: int,
    digestion: Digestion,
    normalize_i_to_l: bool,
    reject_shared_foreign: bool,
) -> ForeignSpeciesEntrapmentGeneration:
    """Construct foreign-species entrapment behavior."""
    return ForeignSpeciesEntrapmentGeneration(
        fold,
        seed,
        digestion,
        normalize_i_to_l,
        reject_shared_foreign,
    )


def _annotation(
    strategy: EntrapmentStrategy,
    fold: int,
    seed: int,
    digestion: Digestion,
    normalize_i_to_l: bool,
    batch: EntrapmentBatch | None,
) -> str:
    note = (
        f"entrapment {strategy.value} fold {fold} seed {seed} with "
        f"protein_fasta {_IMPLEMENTATION_VERSION} "
        f"({digestion.cleavage.pattern}, missed cleavages {digestion.missed_cleavages}, "
        f"length {digestion.min_length}-{digestion.max_length}"
    )
    if normalize_i_to_l:
        note += ", isoleucine normalized to leucine"
    if batch is not None and not batch.is_complete:
        note += (
            f", {batch.complete_proteins} of {batch.source_proteins} proteins at full fold"
            f", {batch.failures} peptides unresolved"
        )
    return f"{note})"


def _digestion_parameters(
    digestion: Digestion,
    *,
    normalize_i_to_l: bool,
) -> dict[str, object]:
    return {
        "enzyme": digestion.cleavage.pattern,
        "missed_cleavages": digestion.missed_cleavages,
        "min_length": digestion.min_length,
        "max_length": digestion.max_length,
        "normalize_i_to_l": normalize_i_to_l,
        "cleavage_residues": ["K", "R"],
    }


def _entrapment_id(identifier: str, fold_index: int, fold: int) -> str:
    if fold == 1:
        return f"{identifier}{_ENTRAPMENT_SUFFIX}"
    return f"{identifier}_{fold_index}{_ENTRAPMENT_SUFFIX}"


def _require_unique_ids(entries: tuple[Entry, ...]) -> None:
    identifiers = [_identifier(description) for description, _ in entries]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Target and foreign FASTA identifiers must be globally unique.")


def _identifier(description: str) -> str:
    return description.partition(" ")[0]
