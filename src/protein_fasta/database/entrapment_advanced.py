"""Optional implementations for biological entrapment generation."""

from __future__ import annotations

import importlib.metadata
from dataclasses import asdict, dataclass
from typing import cast

from fdr_benchmark.entrapment.foreign_species import generate_foreign_species_entrapment
from fdr_benchmark.entrapment.shuffled import generate_shuffled_entrapment
from fdr_benchmark.models import (
    DigestionSpec,
    EntrapmentGenerationResult,
    EntrapmentRequest,
    ExhaustionPolicy,
    FastaRecord,
)

from protein_fasta.database.entrapment import (
    EntrapmentBatch,
    EntrapmentPeptidePair,
    EntrapmentStrategy,
    Entry,
)

ENTRAPMENT_SUFFIX = "_p_target"


@dataclass(frozen=True, slots=True)
class ShuffledEntrapmentGeneration:
    """Generate peptide-shuffled entrapment records."""

    _request: EntrapmentRequest
    _implementation_version: str

    @property
    def strategy(self) -> EntrapmentStrategy:
        """Return the shuffled strategy identity."""
        return EntrapmentStrategy.SHUFFLED

    @property
    def seed(self) -> int:
        """Return the configured seed."""
        return self._request.seed

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]:
        """Apply the configured I/L normalization to every biological entry."""
        return _normalize_entries(entries, self._request.digestion.normalize_i_to_l)

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch:
        """Generate shuffled entries; foreign sources are not used."""
        del foreign_entries
        records = tuple(
            _benchmark_record(description, sequence) for description, sequence in entries
        )
        result = generate_shuffled_entrapment(records, self._request)
        return _batch(result, self.parameters(), len(records))

    def parameters(self) -> dict[str, object]:
        """Return shuffled-entrapment provenance."""
        return _parameters(self.strategy, self._request, self._implementation_version)

    def annotation(self, batch: EntrapmentBatch | None = None) -> str:
        """Render shuffled-entrapment sentinel provenance."""
        return _annotation(self.strategy, self._request, self._implementation_version, batch)


@dataclass(frozen=True, slots=True)
class ForeignSpeciesEntrapmentGeneration:
    """Select entrapment records from a supplied foreign-species database."""

    _request: EntrapmentRequest
    _implementation_version: str

    @property
    def strategy(self) -> EntrapmentStrategy:
        """Return the foreign-species strategy identity."""
        return EntrapmentStrategy.FOREIGN_SPECIES

    @property
    def seed(self) -> int:
        """Return the configured seed."""
        return self._request.seed

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]:
        """Apply the configured I/L normalization to every biological entry."""
        return _normalize_entries(entries, self._request.digestion.normalize_i_to_l)

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch:
        """Select entrapment proteins from the supplied foreign entries."""
        if not foreign_entries:
            raise ValueError(
                "Foreign-species entrapment needs a source database to draw proteins from"
            )
        records = tuple(
            _benchmark_record(description, sequence) for description, sequence in entries
        )
        foreign = tuple(
            _benchmark_record(description, sequence) for description, sequence in foreign_entries
        )
        result = generate_foreign_species_entrapment(records, foreign, self._request)
        return _batch(result, self.parameters(), len(records))

    def parameters(self) -> dict[str, object]:
        """Return foreign-species entrapment provenance."""
        return _parameters(self.strategy, self._request, self._implementation_version)

    def annotation(self, batch: EntrapmentBatch | None = None) -> str:
        """Render foreign-species sentinel provenance."""
        return _annotation(self.strategy, self._request, self._implementation_version, batch)


def make_shuffled_entrapment_generation(
    *,
    fold: int,
    seed: int,
    enzyme: str,
    missed_cleavages: int,
    minimum_length: int,
    maximum_length: int,
    normalize_i_to_l: bool,
    fix_peptide_n_term: bool,
    fix_peptide_c_term: bool,
) -> ShuffledEntrapmentGeneration:
    """Construct peptide-shuffled entrapment behavior from resolved values."""
    request = _request(
        fold=fold,
        seed=seed,
        enzyme=enzyme,
        missed_cleavages=missed_cleavages,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        normalize_i_to_l=normalize_i_to_l,
        fix_peptide_n_term=fix_peptide_n_term,
        fix_peptide_c_term=fix_peptide_c_term,
        reject_shared_foreign=False,
    )
    return ShuffledEntrapmentGeneration(
        request,
        importlib.metadata.version("fdr_benchmark"),
    )


def make_foreign_species_entrapment_generation(
    *,
    fold: int,
    seed: int,
    enzyme: str,
    missed_cleavages: int,
    minimum_length: int,
    maximum_length: int,
    normalize_i_to_l: bool,
    reject_shared_foreign: bool,
) -> ForeignSpeciesEntrapmentGeneration:
    """Construct foreign-species entrapment behavior from resolved values."""
    request = _request(
        fold=fold,
        seed=seed,
        enzyme=enzyme,
        missed_cleavages=missed_cleavages,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        normalize_i_to_l=normalize_i_to_l,
        fix_peptide_n_term=True,
        fix_peptide_c_term=True,
        reject_shared_foreign=reject_shared_foreign,
    )
    return ForeignSpeciesEntrapmentGeneration(
        request,
        importlib.metadata.version("fdr_benchmark"),
    )


def _request(
    *,
    fold: int,
    seed: int,
    enzyme: str,
    missed_cleavages: int,
    minimum_length: int,
    maximum_length: int,
    normalize_i_to_l: bool,
    fix_peptide_n_term: bool,
    fix_peptide_c_term: bool,
    reject_shared_foreign: bool,
) -> EntrapmentRequest:
    return EntrapmentRequest(
        digestion=DigestionSpec(
            enzyme=enzyme,
            missed_cleavages=missed_cleavages,
            min_length=minimum_length,
            max_length=maximum_length,
            normalize_i_to_l=normalize_i_to_l,
        ),
        fold=fold,
        seed=seed,
        suffix=ENTRAPMENT_SUFFIX,
        fix_peptide_n_term=fix_peptide_n_term,
        fix_peptide_c_term=fix_peptide_c_term,
        exhaustion_policy=ExhaustionPolicy.EMIT_PARTIAL,
        reject_shared_foreign_proteins=reject_shared_foreign,
    )


def _batch(
    result: EntrapmentGenerationResult,
    parameters: dict[str, object],
    source_proteins: int,
) -> EntrapmentBatch:
    generated = tuple((record.header, record.sequence) for record in result.records)
    return EntrapmentBatch(
        generated,
        tuple(
            EntrapmentPeptidePair(
                source_id=pair.source_id,
                target_peptide=pair.target_peptide,
                generated_peptide=pair.generated_peptide,
                fold_index=pair.fold_index,
            )
            for pair in result.peptide_pairs
        ),
        parameters,
        result.stats.requested_fold,
        result.stats.achieved_fold,
        result.stats.failures,
        len({failure.source_id for failure in result.failures}),
        source_proteins,
    )


def _normalize_entries(
    entries: tuple[Entry, ...],
    normalize_i_to_l: bool,
) -> tuple[Entry, ...]:
    if not normalize_i_to_l:
        return entries
    return tuple((description, sequence.replace("I", "L")) for description, sequence in entries)


def _benchmark_record(description: str, sequence: str) -> FastaRecord:
    identifier, separator, remainder = description.partition(" ")
    return FastaRecord(identifier, sequence, remainder if separator else "")


def _parameters(
    strategy: EntrapmentStrategy,
    request: EntrapmentRequest,
    implementation_version: str,
) -> dict[str, object]:
    return {
        "strategy": strategy.value,
        "fold": request.fold,
        "seed": request.seed,
        "suffix": ENTRAPMENT_SUFFIX,
        "digestion": cast("dict[str, object]", asdict(request.digestion)),
        "fix_peptide_n_term": request.fix_peptide_n_term,
        "fix_peptide_c_term": request.fix_peptide_c_term,
        "reject_shared_foreign_proteins": request.reject_shared_foreign_proteins,
        "exhaustion_policy": request.exhaustion_policy.value,
        "implementation": "fdr_benchmark",
        "implementation_version": implementation_version,
    }


def _annotation(
    strategy: EntrapmentStrategy,
    request: EntrapmentRequest,
    implementation_version: str,
    batch: EntrapmentBatch | None,
) -> str:
    parameters = _parameters(strategy, request, implementation_version)
    digestion = cast("dict[str, object]", parameters["digestion"])
    note = (
        f"entrapment {strategy.value} fold {request.fold} seed {request.seed} with "
        f"{parameters['implementation']} {parameters['implementation_version']} "
        f"({digestion['enzyme']}, missed cleavages {digestion['missed_cleavages']}, "
        f"length {digestion['min_length']}-{digestion['max_length']}"
    )
    if request.digestion.normalize_i_to_l:
        note += ", isoleucine normalized to leucine"
    if batch is not None and not batch.is_complete:
        note += (
            f", {batch.complete_proteins} of {batch.source_proteins} proteins at full fold"
            f", {batch.failures} peptides unresolved"
        )
    return f"{note})"
