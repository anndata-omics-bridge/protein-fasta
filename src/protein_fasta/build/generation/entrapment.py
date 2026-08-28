"""Generate entrapment FASTA entries.

An entrapment record is a known-false protein placed in the *target* space, so
the rate at which a search identifies one measures the true false discovery
proportion. That is the opposite role from a decoy, which the engine knows about
and estimates FDR from, and the two compose: entrapment records join the target
space before decoys are generated, so they get decoy counterparts too.

Both strategies delegate to the standalone ``fdr_benchmark`` package, exactly as
the shuffle and DecoyPYrat decoy modes do.
"""

from __future__ import annotations

import importlib.metadata
from dataclasses import asdict, dataclass
from typing import Any, Protocol, cast

from fdr_benchmark.entrapment.foreign_species import generate_foreign_species_entrapment
from fdr_benchmark.entrapment.shuffled import generate_shuffled_entrapment
from fdr_benchmark.models import (
    DigestionSpec,
    EntrapmentGenerationResult,
    EntrapmentRequest,
    ExhaustionPolicy,
    FastaRecord,
    PeptidePairRecord,
)

from protein_fasta.analytics_compile import make_digestion
from protein_fasta.schema.build import EntrapmentDocument, EntrapmentStrategy

Entry = tuple[str, str]

ENTRAPMENT_SUFFIX = "_p_target"


# FDRBench digests with one missed cleavage by default, and the collision
# universe is the only thing missed cleavages reach here: the library reconstructs
# proteins from the zero-missed-cleavage partition whatever the request says. One
# is therefore the stricter setting, because an entrapment peptide then also has
# to differ from every target peptide that spans a cleavage site.
DEFAULT_ENTRAPMENT_MISSED_CLEAVAGES = 1


DEFAULT_ENTRAPMENT_SPEC = EntrapmentDocument()


@dataclass(frozen=True, slots=True)
class EntrapmentBatch:
    """Generated entries, the peptide mapping, and algorithm evidence."""

    entries: tuple[Entry, ...]
    peptide_pairs: tuple[PeptidePairRecord, ...]
    parameters: dict[str, Any]
    requested_fold: int
    achieved_fold: int
    failures: int
    proteins_affected: int
    source_proteins: int

    @property
    def is_complete(self) -> bool:
        """Whether every source protein reached the requested fold."""
        return not self.failures and self.achieved_fold >= self.requested_fold

    @property
    def complete_proteins(self) -> int:
        """Source proteins that reached the requested fold in full."""
        return self.source_proteins - self.proteins_affected


class EntrapmentGeneration(Protocol):
    """Generate and describe one compiled entrapment algorithm."""

    @property
    def strategy(self) -> EntrapmentStrategy: ...

    @property
    def seed(self) -> int: ...

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]: ...

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch: ...

    def parameters(self) -> dict[str, Any]: ...

    def annotation(self, batch: EntrapmentBatch | None = None) -> str: ...


@dataclass(frozen=True, slots=True)
class ShuffledEntrapmentGeneration:
    """Generate peptide-shuffled entrapment records."""

    _request: EntrapmentRequest
    _implementation_version: str

    @property
    def strategy(self) -> EntrapmentStrategy:
        return EntrapmentStrategy.SHUFFLED

    @property
    def seed(self) -> int:
        return self._request.seed

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]:
        return _normalize_entries(entries, self._request.digestion.normalize_i_to_l)

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch:
        del foreign_entries
        records = tuple(
            _benchmark_record(description, sequence) for description, sequence in entries
        )
        result = generate_shuffled_entrapment(records, self._request)
        return _batch(result, self.parameters(), len(records))

    def parameters(self) -> dict[str, Any]:
        return _parameters(self.strategy, self._request, self._implementation_version)

    def annotation(self, batch: EntrapmentBatch | None = None) -> str:
        return _annotation(self.strategy, self._request, self._implementation_version, batch)


@dataclass(frozen=True, slots=True)
class ForeignSpeciesEntrapmentGeneration:
    """Select entrapment records from a supplied foreign-species database."""

    _request: EntrapmentRequest
    _implementation_version: str

    @property
    def strategy(self) -> EntrapmentStrategy:
        return EntrapmentStrategy.FOREIGN_SPECIES

    @property
    def seed(self) -> int:
        return self._request.seed

    def normalize(self, entries: tuple[Entry, ...]) -> tuple[Entry, ...]:
        return _normalize_entries(entries, self._request.digestion.normalize_i_to_l)

    def generate(
        self,
        entries: tuple[Entry, ...],
        *,
        foreign_entries: tuple[Entry, ...] = (),
    ) -> EntrapmentBatch:
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

    def parameters(self) -> dict[str, Any]:
        return _parameters(self.strategy, self._request, self._implementation_version)

    def annotation(self, batch: EntrapmentBatch | None = None) -> str:
        return _annotation(self.strategy, self._request, self._implementation_version, batch)


def make_entrapment_generation(spec: EntrapmentDocument, /) -> EntrapmentGeneration:
    """Compile one passive document into behavior selected exactly once."""
    request = _entrapment_request(spec)
    implementation_version = importlib.metadata.version("fdr_benchmark")
    if spec.strategy == EntrapmentStrategy.SHUFFLED:
        return ShuffledEntrapmentGeneration(request, implementation_version)
    if spec.strategy == EntrapmentStrategy.FOREIGN_SPECIES:
        return ForeignSpeciesEntrapmentGeneration(request, implementation_version)
    raise ValueError(f"unsupported entrapment strategy: {spec.strategy!r}")


def generate_entrapment(
    entries: tuple[Entry, ...],
    *,
    spec: EntrapmentDocument = DEFAULT_ENTRAPMENT_SPEC,
    foreign_entries: tuple[Entry, ...] = (),
) -> EntrapmentBatch:
    """Compile and generate entrapment records for compatibility callers."""
    return make_entrapment_generation(spec).generate(
        entries,
        foreign_entries=foreign_entries,
    )


def _batch(
    result: EntrapmentGenerationResult,
    parameters: dict[str, Any],
    source_proteins: int,
) -> EntrapmentBatch:
    """Render generated records as FASTA entries and count what fell short."""
    generated = tuple((record.header, record.sequence) for record in result.records)
    return EntrapmentBatch(
        generated,
        result.peptide_pairs,
        parameters,
        result.stats.requested_fold,
        result.stats.achieved_fold,
        result.stats.failures,
        len({failure.source_id for failure in result.failures}),
        source_proteins,
    )


def _entrapment_request(spec: EntrapmentDocument) -> EntrapmentRequest:
    """Describe this application's digestion to the entrapment generators.

    The collision universe has to be the peptide set a search of this database
    would see, so it follows the configured digestion rather than the generator's
    Pyteomics defaults, as ``_decoypyrat_request`` already does. Missed cleavages
    reach only that universe: ``digest_segments`` reconstructs proteins from the
    zero-missed-cleavage partition whatever this request says.

    ``EMIT_PARTIAL`` rather than the library default ``ERROR``: at human scale a
    minority of low-complexity peptides admit no distinct non-target arrangement,
    and refusing the whole database over them is the DecoyPYrat abort fixed on
    2026-08-11. What was achieved is reported instead, and a build that fell
    short says so in the review and in the sentinel.
    """
    digestion = make_digestion(spec.digestion)
    return EntrapmentRequest(
        digestion=DigestionSpec(
            enzyme=digestion.cleavage.pattern,
            missed_cleavages=spec.digestion.missed_cleavages,
            min_length=spec.digestion.min_length,
            max_length=spec.digestion.max_length,
            normalize_i_to_l=spec.normalize_i_to_l,
        ),
        fold=spec.fold,
        seed=spec.seed,
        suffix=ENTRAPMENT_SUFFIX,
        fix_peptide_n_term=spec.fix_peptide_n_term,
        fix_peptide_c_term=spec.fix_peptide_c_term,
        exhaustion_policy=ExhaustionPolicy.EMIT_PARTIAL,
        reject_shared_foreign_proteins=spec.reject_shared_foreign,
    )


def normalize_entries(entries: tuple[Entry, ...], spec: EntrapmentDocument) -> tuple[Entry, ...]:
    """Apply the spec's sequence normalization to a whole target space.

    FDRBench's ``-I2L`` rewrites the database it emits, not only the entrapment
    half. Doing it to the generated records alone would make them the only
    sequences in the file with no isoleucine.
    """
    return make_entrapment_generation(spec).normalize(entries)


def _normalize_entries(entries: tuple[Entry, ...], normalize_i_to_l: bool) -> tuple[Entry, ...]:
    if not normalize_i_to_l:
        return entries
    return tuple((description, sequence.replace("I", "L")) for description, sequence in entries)


def _benchmark_record(description: str, sequence: str) -> FastaRecord:
    identifier, separator, remainder = description.partition(" ")
    return FastaRecord(identifier, sequence, remainder if separator else "")


def entrapment_parameters(spec: EntrapmentDocument) -> dict[str, Any]:
    """Return machine-readable parameters for sentinel-sidecar provenance."""
    return make_entrapment_generation(spec).parameters()


def _parameters(
    strategy: EntrapmentStrategy,
    request: EntrapmentRequest,
    implementation_version: str,
) -> dict[str, Any]:
    return {
        "strategy": strategy.value,
        "fold": request.fold,
        "seed": request.seed,
        "suffix": ENTRAPMENT_SUFFIX,
        "digestion": asdict(request.digestion),
        "fix_peptide_n_term": request.fix_peptide_n_term,
        "fix_peptide_c_term": request.fix_peptide_c_term,
        "reject_shared_foreign_proteins": request.reject_shared_foreign_proteins,
        "exhaustion_policy": request.exhaustion_policy.value,
        "implementation": "fdr_benchmark",
        "implementation_version": implementation_version,
    }


def entrapment_annotation(spec: EntrapmentDocument, batch: EntrapmentBatch | None = None) -> str:
    """Render compact entrapment provenance for the FASTA sentinel."""
    return make_entrapment_generation(spec).annotation(batch)


def _annotation(
    strategy: EntrapmentStrategy,
    request: EntrapmentRequest,
    implementation_version: str,
    batch: EntrapmentBatch | None,
) -> str:
    parameters = _parameters(strategy, request, implementation_version)
    digestion = cast(dict[str, Any], parameters["digestion"])
    note = (
        f"entrapment {strategy.value} fold {request.fold} seed {request.seed} with "
        f"{parameters['implementation']} {parameters['implementation_version']} "
        f"({digestion['enzyme']}, missed cleavages {digestion['missed_cleavages']}, "
        f"length {digestion['min_length']}-{digestion['max_length']}"
    )
    if request.digestion.normalize_i_to_l:
        # This rewrote every sequence in the file, so it cannot stay implicit.
        note += ", isoleucine normalized to leucine"
    if batch is not None and not batch.is_complete:
        # Not the achieved fold on its own: that is a minimum across proteins, so
        # one hopeless peptide reports zero for a database that is otherwise whole.
        note += (
            f", {batch.complete_proteins} of {batch.source_proteins} proteins at full fold"
            f", {batch.failures} peptides unresolved"
        )
    note += ")"
    return note
