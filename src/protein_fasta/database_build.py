"""Assemble reproducible protein FASTA databases.

The assembly operation is a Python port of prozor's ``create_fgcz_fasta_db``.

Assembly order mirrors prozor:

    [db sentinel] + targets + [contaminants] + [entrapment]

with two FGCZ-convention additions from the newer collection: an ``aa|`` database
sentinel as the first entry, and an ``aa|Cont_...`` section marker before each
contaminant set's block. Entrapment records belong to the biological target
space. Decoy generation is deliberately absent here and consumes the resulting
protein inventory through the separate ``decoy_database`` workflow.

Every input sequence is normalized first -- upper-cased, one trailing stop removed
-- so what is written, what is decoyed, and what deduplication compares are the same
content. Duplicate entries are then removed by id, the first occurrence winning as
in prozor's ``!duplicated(names(...))``, but only when their sequences agree: an id
carrying two different sequences is a disagreement between sources and raises
rather than silently losing one of them.
"""

from __future__ import annotations

import datetime
import importlib.metadata
import os
import tempfile
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from loguru import logger

from protein_fasta.analytics.hashing import FILE_CHECKSUM_VERSION, file_checksum, sequence_hash
from protein_fasta.artifact_io import publish_exclusive, temporary_sibling, write_json_atomic
from protein_fasta.build.metadata import build_section_marker_header, build_sentinel_header
from protein_fasta.build.naming import build_dbname, build_fasta_name
from protein_fasta.database.models import (
    BiologicalDatabase,
    ProteinInventoryEntry,
    SourceRole,
)
from protein_fasta.diagnostics.messages import describe_illegal_residues
from protein_fasta.diagnostics.runtime import DiagnosticRules
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import FastaRecord
from protein_fasta.reading.writer import write_records
from protein_fasta.registry.classification import ContaminantBlockState, classify_record
from protein_fasta.registry.kinds import EntryKind
from protein_fasta.registry.rules import RegistryDiagnosticRules, load_registry_diagnostics
from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.build import (
    BiologicalEntrapmentDocument,
    DatabaseBuildCountsDocument,
    DatabaseBuildEntrapmentEvidenceDocument,
    DatabaseBuildNormalizationDocument,
    DatabaseBuildProfileDocument,
    DatabaseBuildRequestDocument,
    DatabaseBuildResultDocument,
    DatabaseBuildSummaryDocument,
    EffectiveDatabaseBuildDocument,
    EntrapmentDocument,
    EntrapmentStrategy,
    ForeignSpeciesEntrapmentDocument,
    MetadataDocument,
    NamingDocument,
)
from protein_fasta.summary import FastaSummary, summarize_sequences
from protein_fasta.validation.sequence import normalize_sequence

if TYPE_CHECKING:
    import polars as pl

    from protein_fasta.build.generation.entrapment_types import (
        EntrapmentBatch,
        EntrapmentGeneration,
    )

Entry = tuple[str, str]

#: Conflicting ids named in the error before it summarizes the rest. Enough to act
#: on, few enough to read.
_REPORTED_CONFLICTS = 5


def _make_entrapment_generation(spec: EntrapmentDocument) -> EntrapmentGeneration:
    """Load the optional entrapment adapter only when a build selects it."""
    try:
        from protein_fasta.build.generation.entrapment import make_entrapment_generation
    except ModuleNotFoundError as error:
        missing_name = error.name or ""
        if missing_name == "fdr_benchmark" or missing_name.startswith("fdr_benchmark."):
            message = "entrapment generation requires the 'protein-fasta[generation]' extra"
            raise RuntimeError(message) from error
        raise
    return make_entrapment_generation(spec)


@dataclass(frozen=True, slots=True)
class ContaminantBlock:
    """One resolved contaminant module supplied to the build pipeline."""

    name: str
    description: str
    entries: tuple[Entry, ...]


@dataclass(frozen=True, slots=True)
class ProteinSourceProvenance:
    """Prepared-source coordinates retained for one biological identifier."""

    source_order: int
    record_order: int
    source_id: str
    source_role: SourceRole


@dataclass(frozen=True)
class DecoyBuildEvidence:
    """Decoy algorithm identity, parameters, and collision outcomes."""

    mode: str
    seed: int | None
    parameters: dict[str, Any]
    initial_collisions: int
    unresolved_collisions: int
    dropped_peptides: int
    omitted_decoys: int


@dataclass(frozen=True)
class EntrapmentBuildEvidence:
    """Entrapment algorithm identity and achieved multiplicity."""

    strategy: str
    seed: int
    requested_fold: int
    achieved_fold: int
    failures: int
    proteins_affected: int
    source_proteins: int


@dataclass(frozen=True)
class PipelineResult:
    """Outcome of a database build."""

    path: Path
    dbname: str
    n_target: int
    n_contaminant: int
    n_decoy: int
    n_total: int
    contaminant_sets: list[str]
    summary: FastaSummary
    decoy: DecoyBuildEvidence | None = None
    n_entrapment: int = 0
    entrapment: EntrapmentBuildEvidence | None = None
    entrapment_pairs_path: Path | None = None
    # What normalization and deduplication changed on the way in, so a produced
    # file can account for the difference from its sources.
    upper_cased_entries: int = 0
    stop_stripped_entries: int = 0
    duplicates_dropped: int = 0
    database: BiologicalDatabase | None = None


@dataclass(frozen=True, slots=True)
class DatabaseBuildOverrides:
    """Common explicitly supplied values that take precedence over JSON."""

    date: datetime.date | None = None


DEFAULT_BUILD_OVERRIDES = DatabaseBuildOverrides()


@dataclass(frozen=True, slots=True)
class DatabaseBuildExecution:
    """Runtime result and its durable JSON evidence paths."""

    result: PipelineResult
    database: BiologicalDatabase
    document: DatabaseBuildResultDocument
    protein_input_path: Path
    inventory_path: Path
    effective_request_path: Path
    result_path: Path


def _entrapment_build_evidence(
    generation: EntrapmentGeneration | None,
    batch: EntrapmentBatch | None,
) -> EntrapmentBuildEvidence | None:
    if generation is None or batch is None:
        return None
    return EntrapmentBuildEvidence(
        strategy=generation.strategy.value,
        seed=generation.seed,
        requested_fold=batch.requested_fold,
        achieved_fold=batch.achieved_fold,
        failures=batch.failures,
        proteins_affected=batch.proteins_affected,
        source_proteins=batch.source_proteins,
    )


def resolve_database_build(
    profile: DatabaseBuildProfileDocument,
    request: DatabaseBuildRequestDocument,
    *,
    profile_base: Path,
    request_base: Path,
    overrides: DatabaseBuildOverrides = DEFAULT_BUILD_OVERRIDES,
) -> EffectiveDatabaseBuildDocument:
    """Resolve profile, request, and explicit overrides into replayable paths."""
    naming = request.naming or profile.naming
    metadata = request.metadata or profile.metadata
    diagnostics = _resolved_profile_or_request_path(
        profile.diagnostics,
        request.diagnostics,
        request_field_is_set="diagnostics" in request.model_fields_set,
        profile_base=profile_base,
        request_base=request_base,
    )
    return EffectiveDatabaseBuildDocument(
        output_dir=_resolved_path(request.output_dir, request_base),
        date=overrides.date or request.date,
        name_fields=request.name_fields,
        template=request.template or naming.default_dbname,
        naming=naming,
        metadata=metadata,
        diagnostics=diagnostics,
        entrapment=request.entrapment,
        annotation=request.annotation,
        installer=request.installer,
    )


def run_database_build(
    protein_input_path: Path,
    effective: EffectiveDatabaseBuildDocument,
    /,
) -> DatabaseBuildExecution:
    """Read one effective request, build its FASTA, and persist typed evidence."""
    effective.output_dir.mkdir(parents=True, exist_ok=True)
    expected_name = build_fasta_name(
        config=effective.naming,
        template=effective.template,
        date=effective.date,
        decoy=False,
        entrapment=effective.entrapment is not None,
        **effective.name_fields,
    )
    expected_path = effective.output_dir / expected_name
    effective_path = expected_path.with_suffix(f"{expected_path.suffix}.effective.json")
    inventory_path = expected_path.with_suffix(f"{expected_path.suffix}.protein-inventory.parquet")
    pairs_path = expected_path.with_suffix(f"{expected_path.suffix}.entrapment-pairs.parquet")
    result_path = expected_path.with_suffix(f"{expected_path.suffix}.result.json")
    _refuse_existing_build_artifacts(
        expected_path,
        inventory_path,
        result_path,
        *((pairs_path,) if effective.entrapment is not None else ()),
    )
    write_json_atomic(
        effective_path,
        effective.model_dump(mode="json"),
        replace_existing=True,
    )

    input_path = protein_input_path.resolve()
    target_entries, blocks, foreign_entries, source_provenance = _sources_from_protein_input(
        input_path
    )

    diagnostics = load_registry_diagnostics(effective.diagnostics)
    with tempfile.TemporaryDirectory(
        prefix=f".{expected_name}.",
        dir=effective.output_dir,
    ) as staging_name:
        staging_dir = Path(staging_name)
        staged_result = build_database(
            targets=target_entries,
            name_fields=effective.name_fields,
            naming=effective.naming,
            metadata=effective.metadata,
            diagnostics=diagnostics,
            output_dir=staging_dir,
            date=effective.date,
            template=effective.template,
            contaminant_blocks=blocks,
            entrapment_spec=_legacy_entrapment_spec(effective.entrapment),
            foreign_entries=foreign_entries,
            annotation=effective.annotation,
            installer=effective.installer,
            source_provenance=source_provenance,
        )
        if staged_result.database is None:
            raise AssertionError("biological assembly did not return its canonical database")
        database = staged_result.database
        staged_inventory = staging_dir / inventory_path.name
        inventory_rows = write_protein_inventory(staged_result, staged_inventory)
        final_pairs = pairs_path if staged_result.entrapment_pairs_path is not None else None
        result = replace(
            staged_result,
            path=expected_path,
            entrapment_pairs_path=final_pairs,
        )
        result_document = _result_document(
            effective,
            result,
            effective_path,
            input_path,
            inventory_path,
            inventory_rows,
            fasta_artifact_path=staged_result.path,
            inventory_artifact_path=staged_inventory,
            pairs_artifact_path=staged_result.entrapment_pairs_path,
        )
        publications = [(staged_result.path, expected_path), (staged_inventory, inventory_path)]
        if staged_result.entrapment_pairs_path is not None:
            publications.append((staged_result.entrapment_pairs_path, pairs_path))
        published: list[Path] = []
        try:
            for staged, destination in publications:
                publish_exclusive(staged, destination)
                published.append(destination)
            write_json_atomic(
                result_path,
                result_document.model_dump(mode="json"),
                replace_existing=False,
            )
        except BaseException:
            for path in published:
                path.unlink(missing_ok=True)
            raise
    return DatabaseBuildExecution(
        result,
        database,
        result_document,
        input_path,
        inventory_path,
        effective_path,
        result_path,
    )


def run_database_build_from_frame(
    frame: pl.DataFrame,
    effective: EffectiveDatabaseBuildDocument,
    /,
) -> DatabaseBuildExecution:
    """Persist one canonical frame for replay, then run the artifact build use case."""
    from protein_fasta.inventory import validate_protein_input_frame

    validate_protein_input_frame(frame)
    expected_name = build_fasta_name(
        config=effective.naming,
        template=effective.template,
        date=effective.date,
        decoy=False,
        entrapment=effective.entrapment is not None,
        **effective.name_fields,
    )
    input_path = (effective.output_dir / expected_name).with_suffix(
        f"{Path(expected_name).suffix}.protein-input.parquet"
    )
    if input_path.exists():
        raise FileExistsError(f"refusing to replace frame build input: {input_path}")
    effective.output_dir.mkdir(parents=True, exist_ok=True)
    with temporary_sibling(input_path) as staged:
        frame.write_parquet(staged)
        publish_exclusive(staged, input_path)
    try:
        return run_database_build(input_path, effective)
    except BaseException:
        input_path.unlink(missing_ok=True)
        raise


def write_protein_inventory(
    result: PipelineResult,
    path: Path,
) -> int:
    """Write the exact canonical tuple used for one biological FASTA."""
    from protein_fasta.inventory import biological_database_frame

    if result.database is None:
        raise ValueError("pipeline result has no canonical biological database")
    with temporary_sibling(path) as staged:
        biological_database_frame(result.database).write_parquet(staged)
        publish_exclusive(staged, path)
    return len(result.database.entries)


def _refuse_existing_build_artifacts(*paths: Path) -> None:
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace existing biological-build artifacts: {names}")


def _resolved_path(path: Path, base: Path) -> Path:
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _resolved_profile_or_request_path(
    profile_path: Path | None,
    request_path: Path | None,
    *,
    request_field_is_set: bool,
    profile_base: Path,
    request_base: Path,
) -> Path | None:
    if request_field_is_set:
        return None if request_path is None else _resolved_path(request_path, request_base)
    return None if profile_path is None else _resolved_path(profile_path, profile_base)


def _sources_from_protein_input(
    path: Path,
) -> tuple[
    list[Entry],
    tuple[ContaminantBlock, ...],
    list[Entry],
    dict[str, ProteinSourceProvenance],
]:
    """Project one validated canonical input frame into biological source groups."""
    from protein_fasta.inventory import read_protein_input

    frame = read_protein_input(path)
    targets: list[Entry] = []
    foreign: list[Entry] = []
    block_order: list[str] = []
    block_descriptions: dict[str, str] = {}
    block_entries: dict[str, list[Entry]] = {}
    provenance: dict[str, ProteinSourceProvenance] = {}
    for row in frame.iter_rows(named=True):
        entry = (str(row["raw_header"]), str(row["sequence"]))
        role = str(row["role"])
        identifier = str(row["id"])
        provenance.setdefault(
            identifier,
            ProteinSourceProvenance(
                source_order=int(str(row["source_order"])),
                record_order=int(str(row["record_order"])),
                source_id=str(row["source_id"]),
                source_role=cast("SourceRole", role),
            ),
        )
        if role == "target":
            targets.append(entry)
        elif role == "foreign":
            foreign.append(entry)
        elif role == "contaminant":
            block_name = row["block_name"]
            if not isinstance(block_name, str) or not block_name:
                raise ValueError(f"protein-input {path} has a contaminant row without a block name")
            if block_name not in block_entries:
                block_order.append(block_name)
                block_entries[block_name] = []
                description = row["block_description"]
                block_descriptions[block_name] = description if isinstance(description, str) else ""
            block_entries[block_name].append(entry)
        else:
            raise ValueError(f"protein-input {path} has unsupported biological role {role!r}")
    if not targets:
        raise ValueError(f"protein-input {path} contains no target rows")
    blocks = tuple(
        ContaminantBlock(
            name=name,
            description=block_descriptions[name],
            entries=tuple(block_entries[name]),
        )
        for name in block_order
    )
    return targets, blocks, foreign, provenance


def _legacy_entrapment_spec(
    spec: BiologicalEntrapmentDocument | None,
) -> EntrapmentDocument | None:
    """Compile strategy-specific biological storage into the current runtime request."""
    if spec is None:
        return None
    if isinstance(spec, ForeignSpeciesEntrapmentDocument):
        return EntrapmentDocument(
            strategy=EntrapmentStrategy.FOREIGN_SPECIES,
            fold=spec.fold,
            seed=spec.seed,
            digestion=spec.digestion,
            normalize_i_to_l=spec.normalize_i_to_l,
            reject_shared_foreign=spec.reject_shared_foreign,
        )
    return EntrapmentDocument(
        strategy=EntrapmentStrategy.SHUFFLED,
        fold=spec.fold,
        seed=spec.seed,
        digestion=spec.digestion,
        fix_peptide_n_term=spec.fix_peptide_n_term,
        fix_peptide_c_term=spec.fix_peptide_c_term,
        normalize_i_to_l=spec.normalize_i_to_l,
    )


def _artifact(
    path: Path,
    *,
    relative_to: Path,
    schema_name: str,
    schema_version: str,
    row_count: int | None = None,
    recorded_path: Path | None = None,
) -> ArtifactDocument:
    return ArtifactDocument(
        schema_name=schema_name,
        schema_version=schema_version,
        path=Path(os.path.relpath(recorded_path or path, start=relative_to)),
        checksum_version=FILE_CHECKSUM_VERSION,
        checksum=file_checksum(path),
        byte_count=path.stat().st_size,
        row_count=row_count,
    )


def _result_document(
    effective: EffectiveDatabaseBuildDocument,
    result: PipelineResult,
    effective_path: Path,
    input_path: Path,
    inventory_path: Path,
    inventory_rows: int,
    *,
    fasta_artifact_path: Path | None = None,
    inventory_artifact_path: Path | None = None,
    pairs_artifact_path: Path | None = None,
) -> DatabaseBuildResultDocument:
    summary = result.summary
    frequencies = (
        {
            residue: count / summary.total_residues
            for residue, count in summary.aa_frequencies.items()
        }
        if summary.total_residues
        else {}
    )
    sidecars = [
        _artifact(
            effective_path,
            relative_to=effective.output_dir,
            schema_name="effective-database-build",
            schema_version=effective.schema_version,
        )
    ]
    if result.entrapment_pairs_path is not None:
        sidecars.append(
            _artifact(
                pairs_artifact_path or result.entrapment_pairs_path,
                relative_to=effective.output_dir,
                schema_name="entrapment-pairs",
                schema_version="1",
                recorded_path=result.entrapment_pairs_path,
            )
        )
    entrapment = None
    if result.entrapment is not None:
        entrapment = DatabaseBuildEntrapmentEvidenceDocument(
            strategy=result.entrapment.strategy,
            seed=result.entrapment.seed,
            requested_fold=result.entrapment.requested_fold,
            achieved_fold=result.entrapment.achieved_fold,
            failures=result.entrapment.failures,
            proteins_affected=result.entrapment.proteins_affected,
            source_proteins=result.entrapment.source_proteins,
        )
    return DatabaseBuildResultDocument(
        protein_fasta_version=importlib.metadata.version("protein_fasta"),
        effective_request=effective,
        input_artifact=_artifact(
            input_path,
            relative_to=effective.output_dir,
            schema_name="protein-input",
            schema_version="1",
        ),
        biological_fasta=_artifact(
            fasta_artifact_path or result.path,
            relative_to=effective.output_dir,
            schema_name="biological-fasta",
            schema_version="1",
            row_count=result.n_total,
            recorded_path=result.path,
        ),
        protein_inventory=_artifact(
            inventory_artifact_path or inventory_path,
            relative_to=effective.output_dir,
            schema_name="protein-inventory",
            schema_version="1",
            row_count=inventory_rows,
            recorded_path=inventory_path,
        ),
        sidecar_artifacts=tuple(sidecars),
        counts=DatabaseBuildCountsDocument(
            target=result.n_target,
            contaminant=result.n_contaminant,
            entrapment=result.n_entrapment,
            total=result.n_total,
        ),
        normalization=DatabaseBuildNormalizationDocument(
            upper_cased=result.upper_cased_entries,
            terminal_stops_stripped=result.stop_stripped_entries,
            duplicates_dropped=result.duplicates_dropped,
        ),
        summary=DatabaseBuildSummaryDocument(
            n_sequences=summary.n_sequences,
            length_min=summary.length_min,
            length_max=summary.length_max,
            length_mean=summary.length_mean,
            length_q1=summary.length_q1,
            length_median=summary.length_median,
            length_q3=summary.length_q3,
            total_residues=summary.total_residues,
            aa_counts=summary.aa_frequencies,
            aa_frequencies=frequencies,
        ),
        entrapment=entrapment,
    )


def _normalize_entries(
    entries: Iterable[Entry],
    rules: DiagnosticRules,
    *,
    source: str,
) -> tuple[list[Entry], int, int]:
    """Normalize build input, refusing what normalization cannot make valid.

    Returns the normalized entries plus how many were upper-cased and how many lost
    a trailing stop. Runs before deduplication, so two entries that differ only by
    case or a terminator are recognised as the duplicate they are rather than
    reported as a conflict.
    """
    normalized_entries: list[Entry] = []
    upper_cased = 0
    stop_stripped = 0
    for description, sequence in entries:
        normalized = normalize_sequence(sequence)
        upper_cased += normalized.upper_cased
        stop_stripped += normalized.stop_stripped
        entry_id = parse_header(description).id
        if not normalized.sequence:
            raise ValueError(f"{source}: entry {entry_id} has a header but no sequence.")
        illegal = rules.illegal_residues(normalized.sequence)
        if illegal:
            raise ValueError(
                f"{source}: entry {entry_id} has {describe_illegal_residues(illegal)}."
            )
        normalized_entries.append((description, normalized.sequence))
    return normalized_entries, upper_cased, stop_stripped


class ConflictingIdError(ValueError):
    """Two entries share an id but not a sequence, so one of them is wrong.

    prozor -- and this pipeline until now -- kept whichever came first and dropped
    the other silently. Dropping half of a genuine disagreement is a decision no
    build should make on an operator's behalf.
    """

    def __init__(self, entry_ids: list[str]) -> None:
        self.entry_ids = entry_ids
        shown = ", ".join(entry_ids[:_REPORTED_CONFLICTS])
        remainder = len(entry_ids) - _REPORTED_CONFLICTS
        if remainder > 0:
            shown += f", and {remainder:,} more"
        super().__init__(
            f"{len(entry_ids):,} identifier(s) appear more than once with different sequences: {shown}. "
            "Resolve the sources before building."
        )


def _deduplicate_by_id(
    groups: Iterable[Iterable[Entry]],
) -> tuple[list[Entry], list[int], int]:
    """Keep the first entry per id, refusing ids whose sequences disagree.

    Returns the survivors, the survivor count per input group, and how many
    identical duplicates were dropped. An identical duplicate is prozor's case and
    stays a silent-but-counted drop; a conflicting one raises.
    """
    seen: dict[str, str] = {}
    entries: list[Entry] = []
    survivor_counts: list[int] = []
    duplicates_dropped = 0
    conflicts: list[str] = []
    for group in groups:
        kept = 0
        for description, sequence in group:
            entry_id = parse_header(description).id
            if entry_id in seen:
                if seen[entry_id] != sequence:
                    conflicts.append(entry_id)
                else:
                    duplicates_dropped += 1
                continue
            seen[entry_id] = sequence
            entries.append((description, sequence))
            kept += 1
        survivor_counts.append(kept)
    if conflicts:
        raise ConflictingIdError(conflicts)
    return entries, survivor_counts, duplicates_dropped


def build_database(
    *,
    targets: Iterable[Entry],
    name_fields: dict[str, Any],
    naming: NamingDocument,
    metadata: MetadataDocument,
    diagnostics: RegistryDiagnosticRules,
    output_dir: Path,
    date: datetime.date,
    template: str | None = None,
    contaminant_blocks: Iterable[ContaminantBlock] = (),
    entrapment_spec: EntrapmentDocument | None = None,
    foreign_entries: Iterable[Entry] = (),
    annotation: str = "",
    installer: str | None = None,
    source_provenance: Mapping[str, ProteinSourceProvenance] | None = None,
) -> PipelineResult:
    """Build a decoy-free biological FASTA from resolved sources and policy.

    ``name_fields`` are substituted into the chosen naming ``template`` (e.g.
    ``{"project": 42261, "dbn": 1, "description": "human"}``). ``targets`` is
    consumed once; pass a list if you need to reuse it.
    """
    resolved_blocks = list(contaminant_blocks)
    entrapment_generation: EntrapmentGeneration | None = (
        None if entrapment_spec is None else _make_entrapment_generation(entrapment_spec)
    )
    sentinel_cfg = metadata
    if installer is not None:
        sentinel_cfg = sentinel_cfg.model_copy(update={"installer": installer})

    filename = build_fasta_name(
        config=naming,
        template=template,
        date=date,
        decoy=False,
        entrapment=entrapment_spec is not None,
        **name_fields,
    )
    # Sentinel/database identity is the name without the date or decoy marker,
    # matching the collection (e.g. file p42261_db1_..._20260605.fasta -> aa|p42261_db1_...).
    dbname = build_dbname(config=naming, template=template, **name_fields)

    # Normalize every input before anything reads a sequence: what is written, what
    # is decoyed, and what deduplication compares are then all the same content.
    target_list, upper_cased, stop_stripped = _normalize_entries(
        targets,
        diagnostics.rules,
        source="Selected targets",
    )
    target_list = _without_existing_packaging(target_list, diagnostics)
    normalized_blocks: list[ContaminantBlock] = []
    for block in resolved_blocks:
        block_entries, block_upper, block_stops = _normalize_entries(
            block.entries,
            diagnostics.rules,
            source=f"Contaminant set {block.name}",
        )
        upper_cased += block_upper
        stop_stripped += block_stops
        normalized_blocks.append(
            replace(block, entries=tuple(_without_existing_packaging(block_entries, diagnostics)))
        )
    resolved_blocks = normalized_blocks
    normalized_foreign, foreign_upper, foreign_stops = _normalize_entries(
        foreign_entries,
        diagnostics.rules,
        source="Foreign entrapment source",
    )
    upper_cased += foreign_upper
    stop_stripped += foreign_stops
    foreign_entries = normalized_foreign
    foreign_entries = _without_existing_packaging(foreign_entries, diagnostics)
    if entrapment_generation is not None:
        target_list = list(entrapment_generation.normalize(tuple(target_list)))
        resolved_blocks = [
            replace(block, entries=entrapment_generation.normalize(block.entries))
            for block in resolved_blocks
        ]
        foreign_entries = entrapment_generation.normalize(tuple(foreign_entries))
    contaminant_entries = [entry for block in resolved_blocks for entry in block.entries]

    biological_source, _, _ = _deduplicate_by_id([target_list, contaminant_entries])
    entrapment_batch: EntrapmentBatch | None = None
    entrapment_entries: tuple[Entry, ...] = ()
    if entrapment_generation is not None:
        entrapment_batch = entrapment_generation.generate(
            tuple(biological_source),
            foreign_entries=tuple(foreign_entries),
        )
        entrapment_entries = entrapment_batch.entries
        entrapment_note = entrapment_generation.annotation(entrapment_batch)
        annotation = f"{annotation}; {entrapment_note}" if annotation else entrapment_note

    # Assemble in prozor order with marker-delimited contaminant blocks, deduplicating
    # by id across the whole file (first occurrence wins, matching prozor's
    # !duplicated). The category counts are the survivors of dedup, so they always
    # reconcile with the written file — important for merged proteomes that share
    # accessions (the multi-select feature's core case).
    sentinel_header = build_sentinel_header(dbname, annotation, date, sentinel_cfg)
    groups: list[Iterable[Entry]] = [
        [(sentinel_header, sentinel_cfg.body_sequence)],
        target_list,
    ]
    for block in resolved_blocks:
        marker_header = build_section_marker_header(
            f"Cont_{block.name}",
            block.description,
            sentinel_cfg,
        )
        groups.append([(marker_header, sentinel_cfg.marker_body_sequence)])
        groups.append(block.entries)
    groups.append(entrapment_entries)
    final, _, duplicates_dropped = _deduplicate_by_id(groups)
    database = BiologicalDatabase(
        _protein_inventory_entries(
            final,
            diagnostics,
            source_provenance or {},
            None if entrapment_generation is None else entrapment_generation.strategy.value,
        )
    )
    n_target = sum(entry.kind == EntryKind.TARGET.value for entry in database.entries)
    n_contaminant = sum(entry.kind == EntryKind.CONTAMINANT.value for entry in database.entries)
    n_entrapment = sum(entry.kind == EntryKind.ENTRAPMENT.value for entry in database.entries)
    n_decoy = 0
    out_path = output_dir / filename
    with temporary_sibling(out_path) as staged_fasta:
        write_records(
            (
                FastaRecord(raw_header=entry.raw_header, sequence=entry.sequence)
                for entry in database.entries
            ),
            staged_fasta,
        )
        publish_exclusive(staged_fasta, out_path)

    pairs_path: Path | None = None
    if entrapment_batch is not None and entrapment_batch.peptide_pairs:
        from protein_fasta.inventory import entrapment_pair_frame

        # Beside the FASTA and named for it: the later FDP evaluation has nothing
        # to join search results on without this mapping, and a sidecar travels
        # with the database the way the build manifest already does.
        #
        # Only when there are pairs. Foreign-species selection produces none,
        # because a foreign protein is not peptide-paired with the target it
        # entraps, and a header-only file would claim a pairing that does not
        # exist -- which a paired FDP estimate would then silently be run against.
        pairs_path = out_path.with_suffix(f"{out_path.suffix}.entrapment-pairs.parquet")
        with temporary_sibling(pairs_path) as staged_pairs:
            entrapment_pair_frame(
                [
                    {
                        "source_id": pair.source_id,
                        "target_peptide": pair.target_peptide,
                        "generated_peptide": pair.generated_peptide,
                        "fold_index": pair.fold_index,
                    }
                    for pair in entrapment_batch.peptide_pairs
                ]
            ).write_parquet(staged_pairs)
            publish_exclusive(staged_pairs, pairs_path)

    summary = summarize_sequences(_scientific_sequences(final, diagnostics))
    logger.info(
        "built {}: {} target + {} contaminant + {} entrapment + {} decoy = {} entries",
        out_path.name,
        n_target,
        n_contaminant,
        n_entrapment,
        n_decoy,
        len(final),
    )
    entrapment_evidence = _entrapment_build_evidence(entrapment_generation, entrapment_batch)
    return PipelineResult(
        path=out_path,
        dbname=dbname,
        n_target=n_target,
        n_contaminant=n_contaminant,
        n_decoy=n_decoy,
        n_total=len(final),
        contaminant_sets=[block.name for block in resolved_blocks],
        summary=summary,
        decoy=None,
        n_entrapment=n_entrapment,
        entrapment=entrapment_evidence,
        entrapment_pairs_path=pairs_path,
        upper_cased_entries=upper_cased,
        stop_stripped_entries=stop_stripped,
        duplicates_dropped=duplicates_dropped,
        database=database,
    )


def _without_existing_packaging(
    records: Iterable[Entry],
    diagnostics: RegistryDiagnosticRules,
) -> list[Entry]:
    """Remove pre-existing sentinel and decoy records from a biological source."""
    biological: list[Entry] = []
    block_state: ContaminantBlockState | None = None
    for raw_header, sequence in records:
        header = parse_header(raw_header)
        _, classifications = diagnostics.rules.diagnose_identifier(header.id)
        kind, _group, block_state = classify_record(
            raw_header,
            classifications,
            block_state,
            diagnostics.decoy_prefix,
        )
        if kind not in {EntryKind.SENTINEL, EntryKind.DECOY}:
            biological.append((raw_header, sequence))
    return biological


def _protein_inventory_entries(
    records: Iterable[Entry],
    diagnostics: RegistryDiagnosticRules,
    provenance_by_id: Mapping[str, ProteinSourceProvenance],
    entrapment_strategy: str | None,
) -> tuple[ProteinInventoryEntry, ...]:
    """Project final biological records once into immutable inventory entries."""
    entries: list[ProteinInventoryEntry] = []
    block_state: ContaminantBlockState | None = None
    for order, (raw_header, sequence) in enumerate(records):
        header = parse_header(raw_header)
        _, classifications = diagnostics.rules.diagnose_identifier(header.id)
        kind, contaminant_group, block_state = classify_record(
            raw_header,
            classifications,
            block_state,
            diagnostics.decoy_prefix,
        )
        if kind is EntryKind.DECOY:
            raise ValueError(f"biological assembly classified {header.id!r} as a decoy")
        provenance = provenance_by_id.get(header.id)
        generated_entrapment = provenance is None and kind is EntryKind.ENTRAPMENT
        entries.append(
            ProteinInventoryEntry(
                final_order=order,
                raw_header=raw_header,
                identifier=header.id,
                description=header.description,
                sequence=sequence,
                kind=kind.value,
                contaminant_group=contaminant_group,
                sequence_hash=sequence_hash(sequence).hex(),
                entrapment_strategy=(entrapment_strategy if kind is EntryKind.ENTRAPMENT else None),
                source_order=None if provenance is None else provenance.source_order,
                record_order=None if provenance is None else provenance.record_order,
                source_id=(
                    "generated-entrapment"
                    if generated_entrapment
                    else None
                    if provenance is None
                    else provenance.source_id
                ),
                source_role=(
                    "entrapment"
                    if generated_entrapment
                    else None
                    if provenance is None
                    else provenance.source_role
                ),
            )
        )
    return tuple(entries)


def _scientific_sequences(
    entries: Iterable[Entry],
    diagnostics: RegistryDiagnosticRules,
) -> Iterable[str]:
    """Yield biological sequences while excluding sentinels and section markers."""
    block_state: ContaminantBlockState | None = None
    for raw_header, sequence in entries:
        header = parse_header(raw_header)
        _, classifications = diagnostics.rules.diagnose_identifier(header.id)
        kind, _group, block_state = classify_record(
            raw_header,
            classifications,
            block_state,
            diagnostics.decoy_prefix,
        )
        if kind in {EntryKind.TARGET, EntryKind.CONTAMINANT, EntryKind.ENTRAPMENT}:
            yield sequence
