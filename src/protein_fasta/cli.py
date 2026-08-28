"""Command-line inspection of protein FASTA databases."""

from __future__ import annotations

import datetime
import sys
from collections.abc import Callable
from pathlib import Path

import polars as pl
from cyclopts import App
from loguru import logger

from protein_fasta.analytics.clustering import ClusteringMetric
from protein_fasta.analytics.digestion import digest_sequence
from protein_fasta.analytics.hashing import (
    CONTENT_FINGERPRINT_VERSION,
    FILE_CHECKSUM_VERSION,
    content_fingerprint,
    file_checksum,
    sequence_hash,
)
from protein_fasta.analytics_compile import make_digestion
from protein_fasta.candidate_analysis import (
    resolve_candidate_request,
    run_candidate_analysis,
)
from protein_fasta.compile import make_diagnostic_rules
from protein_fasta.database_build import (
    DatabaseBuildOverrides,
    resolve_database_build,
    run_database_build,
)
from protein_fasta.decoy_database import (
    DecoyOverrides,
    resolve_decoy_request,
    run_decoy_generation,
)
from protein_fasta.decoy_report import resolve_decoy_report_request, run_decoy_report
from protein_fasta.diagnostic_summary import (
    ProteinDiagnosticsSummary,
    summarize_protein_diagnostics,
)
from protein_fasta.documents import (
    load_builtin_database_build_profile,
    load_builtin_diagnostic_document,
    load_builtin_entry_classifier_document,
    load_builtin_header_format_catalog,
    load_candidate_request,
    load_database_build_profile,
    load_database_build_request,
    load_decoy_report_request,
    load_decoy_request,
    load_derived_protein_input_request,
    load_digestion_document,
    load_entry_classifier_document,
    load_header_format_catalog,
    load_peptide_build_request,
    load_peptide_comparison_request,
    load_protein_input_request,
    load_registry_document,
    load_uniprot_catalog_request,
    load_uniprot_download_request,
)
from protein_fasta.frame import (
    read_basic_protein_frame,
    read_configured_protein_frame,
    read_header_format_diagnostics_frame,
    read_protein_frame,
    read_strict_configured_protein_frame,
    read_strict_protein_frame,
)
from protein_fasta.peptide_workflow import (
    resolve_peptide_build_request,
    resolve_peptide_comparison_request,
    run_peptide_build,
    run_peptide_comparison,
)
from protein_fasta.protein_input import (
    resolve_derived_protein_input_request,
    resolve_protein_input_request,
    run_derived_protein_input_preparation,
    run_protein_input_preparation,
)
from protein_fasta.record import iter_protein_diagnostics, iter_proteins
from protein_fasta.registry.backend import factory as registry_backends
from protein_fasta.registry.backend.base import RegistryConnection
from protein_fasta.registry.clustering import cluster_registered_databases
from protein_fasta.registry.comparisons import compare_database
from protein_fasta.registry.export import query_similarity_data, write_similarity_exports
from protein_fasta.registry.indexing import (
    SCHEMA_VERSION,
    RegistryRecord,
    RegistrySchemaError,
    RejectedFasta,
    connect_registry,
    list_databases,
    rebuild_registry,
    update_registry,
)
from protein_fasta.registry.kinds import EntryKind
from protein_fasta.registry_workflow import index_database_inventory, make_registry_settings
from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.registry import RegistryDocument
from protein_fasta.schema.uniprot import UniProtCatalogRequestDocument
from protein_fasta.uniprot_catalog import (
    filter_uniprot_catalog,
    read_uniprot_catalog,
    sync_uniprot_catalog,
)
from protein_fasta.uniprot_download import (
    UniProtDownloadOverrides,
    resolve_uniprot_download,
    run_uniprot_download,
)

type _TableWriter = Callable[[pl.DataFrame, Path], None]

app = App(
    name="protein-fasta",
    help="Inspect a protein FASTA database as tables or aggregate diagnostics.",
    help_on_error=True,
)


def _write_csv(frame: pl.DataFrame, path: Path) -> None:
    frame.write_csv(path)


def _write_tsv(frame: pl.DataFrame, path: Path) -> None:
    frame.write_csv(path, separator="\t")


def _write_excel(frame: pl.DataFrame, path: Path) -> None:
    frame.write_excel(path)


def _write_parquet(frame: pl.DataFrame, path: Path) -> None:
    frame.write_parquet(path)


_WRITERS: dict[str, _TableWriter] = {
    ".csv": _write_csv,
    ".parquet": _write_parquet,
    ".tsv": _write_tsv,
    ".xlsx": _write_excel,
}


def _writer_for(path: Path) -> _TableWriter:
    writer = _WRITERS.get(path.suffix.lower())
    if writer is None:
        supported = ", ".join(sorted(_WRITERS))
        raise ValueError(f"unsupported table suffix {path.suffix!r}; choose one of: {supported}")
    return writer


def _cli_path(path: Path, /) -> Path:
    """Resolve one command-line path against the caller's working directory."""
    return path.resolve()


def _export_frame(
    frame: pl.DataFrame,
    table_path: Path,
    *,
    exclude_sequence: bool,
) -> None:
    if exclude_sequence:
        frame = frame.drop("sequence")
    _writer_for(table_path)(frame, table_path)
    logger.info(
        "Wrote {} FASTA records and {} columns to {}",
        frame.height,
        frame.width,
        table_path,
    )


def _with_row_checksums(frame: pl.DataFrame) -> pl.DataFrame:
    identifiers: list[str] = frame.get_column("id").to_list()
    sequences: list[str] = frame.get_column("sequence").to_list()
    sequence_hashes = [sequence_hash(sequence) for sequence in sequences]
    checksum_frame = pl.DataFrame(
        {
            "sequence_hash": [hashed_sequence.hex() for hashed_sequence in sequence_hashes],
            "id_sequence_fingerprint": [
                content_fingerprint(((identifier, hashed_sequence),))
                for identifier, hashed_sequence in zip(
                    identifiers,
                    sequence_hashes,
                    strict=True,
                )
            ],
        },
        schema={
            "sequence_hash": pl.String,
            "id_sequence_fingerprint": pl.String,
        },
    )
    return pl.concat((frame, checksum_frame), how="horizontal_extend")


@app.command
def table(
    fasta_path: Path,
    table_path: Path,
    *,
    sequence: bool = True,
    strict: bool = False,
    checksums: bool = False,
) -> None:
    """Export one FASTA database to CSV, TSV, XLSX, or Parquet.

    By default, each row accepted by exactly one built-in parser is enriched and every other row
    retains its base values. Strict mode enriches only when one parser accepts the complete file.
    Use ``--no-sequence`` when the table is intended for inspection rather than sequence work.
    ``--checksums`` adds a normalized sequence hash and an ID/sequence fingerprint per row.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
        table_path: Output path ending in ``.csv``, ``.tsv``, ``.xlsx``, or ``.parquet``.
        sequence: Include the normalized sequence column in the output.
        strict: Require one parser to accept every row before adding enrichment columns.
        checksums: Include per-row sequence hashes and ID/sequence fingerprints.
    """
    reader = read_strict_protein_frame if strict else read_protein_frame
    frame = reader(fasta_path)
    if checksums:
        frame = _with_row_checksums(frame)
    _export_frame(
        frame,
        table_path,
        exclude_sequence=not sequence,
    )


@app.command
def basic(
    fasta_path: Path,
    table_path: Path,
    *,
    sequence: bool = True,
) -> None:
    """Export only the stable id, description, and sequence columns.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
        table_path: Output path ending in ``.csv``, ``.tsv``, ``.xlsx``, or ``.parquet``.
        sequence: Include the normalized sequence column in the output.
    """
    _export_frame(
        read_basic_protein_frame(fasta_path),
        table_path,
        exclude_sequence=not sequence,
    )


@app.command
def configured(
    fasta_path: Path,
    table_path: Path,
    *,
    rules: tuple[Path, ...],
    classifiers: Path,
    sequence: bool = True,
    strict: bool = False,
) -> None:
    """Export a table using explicit JSON parser rules.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
        table_path: Output path ending in ``.csv``, ``.tsv``, ``.xlsx``, or ``.parquet``.
        rules: One or more database-format JSON documents.
        classifiers: Entry-classifier JSON document.
        sequence: Include the normalized sequence column in the output.
        strict: Require one parser to accept every row before adding enrichment columns.
    """
    reader = read_strict_configured_protein_frame if strict else read_configured_protein_frame
    frame = reader(
        fasta_path,
        load_header_format_catalog(rules),
        load_entry_classifier_document(classifiers),
    )
    _export_frame(frame, table_path, exclude_sequence=not sequence)


@app.command
def formats(fasta_path: Path, table_path: Path) -> None:
    """Export one parser-selection diagnostic row per built-in database format.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
        table_path: Output path ending in ``.csv``, ``.tsv``, ``.xlsx``, or ``.parquet``.
    """
    frame = read_header_format_diagnostics_frame(
        fasta_path,
        load_builtin_header_format_catalog(),
        load_builtin_entry_classifier_document(),
    )
    _writer_for(table_path)(frame, table_path)
    logger.info(
        "Wrote {} format diagnostics and {} columns to {}",
        frame.height,
        frame.width,
        table_path,
    )


def _render_counts(counts: dict[str, int]) -> str:
    return ", ".join(f"{name}={count}" for name, count in counts.items()) or "none"


def _render_frequencies(counts: dict[str, int], total: int) -> str:
    if not total:
        return "none"
    return ", ".join(f"{name}={count} ({count / total:.2%})" for name, count in counts.items())


def _render_combinations(counts: dict[tuple[str, ...], int]) -> str:
    return (
        ", ".join(
            f"{'+'.join(names) if names else 'unclassified'}={count}"
            for names, count in counts.items()
        )
        or "none"
    )


def _log_diagnostics_summary(summary: ProteinDiagnosticsSummary) -> None:
    proteins = summary.proteins
    logger.info("Records: {}; total residues: {}", proteins.n_sequences, proteins.total_residues)
    logger.info(
        "Sequence lengths: min={}; q1={}; median={}; mean={:.2f}; q3={}; max={}",
        proteins.length_min,
        proteins.length_q1,
        proteins.length_median,
        proteins.length_mean,
        proteins.length_q3,
        proteins.length_max,
    )
    logger.info(
        "Amino-acid frequencies: {}",
        _render_frequencies(proteins.aa_frequencies, proteins.total_residues),
    )
    logger.info("Identifier namespaces: {}", _render_counts(summary.namespace_counts))
    logger.info("Classifications: {}", _render_counts(summary.classification_counts))
    logger.info(
        "Classification combinations: {}",
        _render_combinations(summary.classification_combination_counts),
    )
    logger.info(
        "Normalization changes: upper-cased={}; terminal-stop-stripped={}",
        summary.upper_cased_count,
        summary.stop_stripped_count,
    )
    logger.info(
        "Illegal sequences: {}; illegal residues: {}",
        summary.illegal_sequence_count,
        _render_counts(summary.illegal_residue_counts),
    )


@app.command
def diagnostics(fasta_path: Path) -> None:
    """Summarize built-in record diagnostics for one FASTA database.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
    """
    rules = make_diagnostic_rules(
        load_builtin_diagnostic_document(),
        load_builtin_entry_classifier_document(),
    )
    summary = summarize_protein_diagnostics(iter_protein_diagnostics(fasta_path, rules))
    _log_diagnostics_summary(summary)


@app.command
def digest(
    fasta_path: Path,
    table_path: Path,
    *,
    config: Path | None = None,
) -> None:
    """Export theoretical peptides and missed-cleavage counts.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
        table_path: Output path ending in ``.csv``, ``.tsv``, ``.xlsx``, or ``.parquet``.
        config: Optional digestion JSON; packaged trypsin defaults are used when omitted.
    """
    document = load_digestion_document(config) if config is not None else DigestionDocument()
    digestion = make_digestion(document)
    rows = [
        {
            "protein_id": protein.id,
            "peptide": peptide.sequence,
            "missed_cleavages": peptide.missed_cleavages,
        }
        for protein in iter_proteins(fasta_path)
        for peptide in digest_sequence(protein.sequence, digestion)
    ]
    frame = pl.DataFrame(
        rows,
        schema={"protein_id": pl.String, "peptide": pl.String, "missed_cleavages": pl.UInt8},
    )
    _writer_for(table_path)(frame, table_path)
    logger.info("Wrote {} theoretical peptides to {}", frame.height, table_path)


@app.command
def peptides(inventory_path: Path, request_path: Path) -> None:
    """Build canonical peptide, mapping, and unique FASTA artifacts.

    Args:
        inventory_path: Canonical biological or search inventory Parquet.
        request_path: Peptide-build request JSON.
    """
    effective = resolve_peptide_build_request(
        load_peptide_build_request(request_path),
        request_base=request_path.resolve().parent,
    )
    execution = run_peptide_build(_cli_path(inventory_path), effective)
    logger.info("Wrote peptide-build evidence to {}", execution.result_path)


@app.command
def pepcompare(
    peptides_a: Path,
    peptides_b: Path,
    request_path: Path,
) -> None:
    """Compare two canonical peptide inventories exactly.

    Args:
        peptides_a: First canonical peptides Parquet.
        peptides_b: Second canonical peptides Parquet.
        request_path: Peptide-comparison request JSON.
    """
    effective = resolve_peptide_comparison_request(
        load_peptide_comparison_request(request_path),
        request_base=request_path.resolve().parent,
    )
    execution = run_peptide_comparison(
        _cli_path(peptides_a),
        _cli_path(peptides_b),
        effective,
    )
    logger.info("Wrote peptide-comparison evidence to {}", execution.result_path)


@app.command
def checksum(fasta_path: Path) -> None:
    """Report exact-file MD5 and normalized protein-content BLAKE2b fingerprints.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
    """
    exact_file = file_checksum(fasta_path)
    normalized_content = content_fingerprint(
        (protein.id, sequence_hash(protein.sequence)) for protein in iter_proteins(fasta_path)
    )
    logger.info("Exact file ({}): {}", FILE_CHECKSUM_VERSION, exact_file)
    logger.info(
        "Normalized ID/sequence content ({}): {}",
        CONTENT_FINGERPRINT_VERSION,
        normalized_content,
    )


@app.command(name="uniprot-catalog")
def uniprot_catalog(
    request_path: Path,
    *,
    output_dir: Path | None = None,
    timeout_seconds: float | None = None,
) -> None:
    """Synchronize one immutable local UniProt proteome-catalog snapshot.

    Args:
        request_path: UniProt catalog request JSON.
        output_dir: Optional explicit output-directory override.
        timeout_seconds: Optional HTTP timeout override.
    """
    request = load_uniprot_catalog_request(request_path)
    updates: dict[str, object] = {}
    if output_dir is not None:
        updates["output_dir"] = _cli_path(output_dir)
    if timeout_seconds is not None:
        updates["timeout_seconds"] = timeout_seconds
    if updates:
        request = UniProtCatalogRequestDocument.model_validate(request.model_dump() | updates)
    execution = sync_uniprot_catalog(
        request,
        request_base=request_path.resolve().parent,
    )
    logger.info(
        "Wrote {} UniProt proteomes to {}",
        execution.frame.height,
        execution.catalog_path,
    )
    logger.info("Wrote catalog evidence to {}", execution.result_path)


@app.command(name="uniprot-proteomes")
def uniprot_proteomes(
    catalog_path: Path,
    table_path: Path,
    *,
    query: str = "",
) -> None:
    """Filter one local UniProt catalog without contacting the provider.

    Args:
        catalog_path: Existing canonical UniProt catalog Parquet.
        table_path: CSV, TSV, XLSX, or Parquet output.
        query: Optional case-insensitive proteome, organism, or taxid filter.
    """
    frame = filter_uniprot_catalog(read_uniprot_catalog(catalog_path), query)
    _writer_for(table_path)(frame, table_path)
    logger.info("Wrote {} UniProt proteomes to {}", frame.height, table_path)


@app.command(name="uniprot-download")
def uniprot_download(
    request_path: Path,
    *,
    output_fasta: Path | None = None,
    timeout_seconds: float | None = None,
) -> None:
    """Acquire one UniProt proteome FASTA without building or adding decoys.

    Args:
        request_path: UniProt download request JSON.
        output_fasta: Optional explicit FASTA destination override.
        timeout_seconds: Optional HTTP timeout override.
    """
    request = load_uniprot_download_request(request_path)
    effective = resolve_uniprot_download(
        request,
        request_base=request_path.resolve().parent,
        overrides=UniProtDownloadOverrides(
            output_fasta=None if output_fasta is None else _cli_path(output_fasta),
            timeout_seconds=timeout_seconds,
        ),
    )
    execution = run_uniprot_download(effective)
    logger.info("Wrote UniProt download evidence to {}", execution.result_path)


@app.command
def prepare(request_path: Path) -> None:
    """Prepare ordered FASTA sources as the canonical protein-input Parquet.

    Args:
        request_path: Protein-input request JSON containing source roles and output path.
    """
    request = load_protein_input_request(request_path)
    effective = resolve_protein_input_request(
        request,
        request_base=request_path.resolve().parent,
    )
    execution = run_protein_input_preparation(effective)
    logger.info("Wrote protein-input evidence to {}", execution.result_path)


@app.command(name="derive-input")
def derive_input(request_path: Path) -> None:
    """Prepare clean source rows from an existing biological or search inventory.

    The source retains target and contaminant proteins. Existing sentinel,
    section-marker, entrapment, and decoy rows are excluded. An optional second
    inventory supplies foreign proteins for a subsequent entrapment build.

    Args:
        request_path: Derived protein-input request JSON with inventory paths and output.
    """
    request = load_derived_protein_input_request(request_path)
    effective = resolve_derived_protein_input_request(
        request,
        request_base=request_path.resolve().parent,
    )
    execution = run_derived_protein_input_preparation(effective)
    logger.info("Wrote derived protein-input evidence to {}", execution.result_path)


@app.command
def decoy(
    biological_inventory: Path,
    request_path: Path,
    *,
    output_fasta: Path | None = None,
    decoy_prefix: str | None = None,
) -> None:
    """Generate one search database from a biological protein inventory.

    Args:
        biological_inventory: Existing decoy-free protein-inventory Parquet.
        request_path: Decoy request JSON containing exactly one strategy.
        output_fasta: Optional explicit search-FASTA destination override.
        decoy_prefix: Optional explicit generated-header prefix override.
    """
    request = load_decoy_request(request_path)
    effective = resolve_decoy_request(
        request,
        request_base=request_path.resolve().parent,
        overrides=DecoyOverrides(
            output_fasta=None if output_fasta is None else _cli_path(output_fasta),
            decoy_prefix=decoy_prefix,
        ),
    )
    execution = run_decoy_generation(_cli_path(biological_inventory), effective)
    logger.info("Wrote decoy-generation evidence to {}", execution.result_path)


@app.command(name="decoy-report")
def decoy_report(biological_inventory: Path, request_path: Path) -> None:
    """Compare requested decoy methods at peptide level.

    Args:
        biological_inventory: Existing decoy-free protein-inventory Parquet.
        request_path: Decoy-method comparison request JSON.
    """
    effective = resolve_decoy_report_request(
        load_decoy_report_request(request_path),
        request_base=request_path.resolve().parent,
    )
    execution = run_decoy_report(_cli_path(biological_inventory), effective)
    logger.info("Wrote decoy-method evidence to {}", execution.result_path)


@app.command
def build(
    protein_input: Path,
    config: Path,
    *,
    profile: Path | None = None,
    date: datetime.date | None = None,
) -> None:
    """Build one reproducible biological FASTA from prepared protein rows.

    Relative request paths resolve beside the request file. Profile defaults are
    overridden by request values and then by explicitly supplied CLI options. The
    effective request is written before sequence work starts; the final result JSON
    records artifacts, checksums, counts, summaries, normalization, and entrapment.
    Decoy generation is the separate ``decoy`` command.

    Args:
        protein_input: Canonical decoy-free protein-input Parquet.
        config: Per-run database-build request JSON.
        profile: Optional portable defaults JSON; packaged FGCZ defaults are used otherwise.
        date: Optional build-date override.
    """
    request = load_database_build_request(config)
    if profile is None:
        build_profile = load_builtin_database_build_profile()
        profile_base = config.parent
    else:
        build_profile = load_database_build_profile(profile)
        profile_base = profile.parent
    effective = resolve_database_build(
        build_profile,
        request,
        profile_base=profile_base,
        request_base=config.parent,
        overrides=DatabaseBuildOverrides(date=date),
    )
    execution = run_database_build(_cli_path(protein_input), effective)
    logger.info("Built {} entries -> {}", execution.result.n_total, execution.result.path)
    logger.info("Effective request -> {}", execution.effective_request_path)
    logger.info("Result -> {}", execution.result_path)


@app.command
def index(
    fasta_dir: Path,
    registry_path: Path,
    *,
    config: Path | None = None,
    recursive: bool = False,
    full: bool = False,
    force: bool = False,
    prune: bool = False,
) -> None:
    """Build or incrementally update a SQLite or DuckDB FASTA registry.

    Args:
        fasta_dir: Directory containing protein FASTA files.
        registry_path: New ``.sqlite3`` or ``.duckdb`` registry file.
        config: Optional registry-policy JSON document.
        recursive: Search subdirectories too.
        full: Rebuild the complete registry atomically.
        force: Reindex unchanged files during an incremental update.
        prune: Remove missing files during an incremental update.
    """
    document = load_registry_document(config) if config is not None else RegistryDocument()
    path_backend = registry_backends.backend_for_path(registry_path)
    if path_backend != document.registry.backend:
        raise ValueError(
            f"registry path selects {path_backend!r}, but configuration selects "
            f"{document.registry.backend!r}"
        )
    rejections: list[RejectedFasta] = []
    settings = make_registry_settings(
        document,
        fasta_root=fasta_dir,
        registry_path=registry_path,
    )
    if full or not registry_path.exists():
        records = rebuild_registry(
            fasta_dir,
            registry_path,
            settings,
            recursive=recursive,
            rejections=rejections,
        )
    else:
        records = update_registry(
            fasta_dir,
            registry_path,
            settings,
            force=force,
            prune=prune,
            recursive=recursive,
            rejections=rejections,
        )
    logger.info(
        "Indexed {} databases into {}; rejected {} files",
        len(records),
        registry_path,
        len(rejections),
    )
    for rejection in rejections:
        logger.warning("Rejected {}: {}", rejection.path, rejection.reason)


@app.command(name="index-inventory")
def index_inventory(
    inventory_path: Path,
    registry_path: Path,
    *,
    config: Path | None = None,
    label: str | None = None,
) -> None:
    """Index one canonical biological or search inventory directly.

    Args:
        inventory_path: Canonical protein- or search-inventory Parquet.
        registry_path: New or existing ``.sqlite3`` or ``.duckdb`` registry.
        config: Optional registry-policy JSON document.
        label: Optional database filename stored in the registry.
    """
    document = load_registry_document(config) if config is not None else RegistryDocument()
    record = index_database_inventory(
        _cli_path(inventory_path),
        _cli_path(registry_path),
        document,
        label=label,
    )
    logger.info("Indexed {} entries from {}", record.entry_count, inventory_path)


@app.command
def candidate(
    inventory_path: Path,
    registry_path: Path,
    request_path: Path,
    *,
    config: Path | None = None,
) -> None:
    """Compare a biological or search inventory without registering it.

    Args:
        inventory_path: Canonical protein- or search-inventory Parquet.
        registry_path: Existing SQLite or DuckDB registry.
        request_path: Candidate-review request JSON.
        config: Optional registry-policy JSON matching the indexed registry.
    """
    request = load_candidate_request(request_path)
    effective = resolve_candidate_request(
        request,
        request_base=request_path.resolve().parent,
    )
    document = load_registry_document(config) if config is not None else RegistryDocument()
    execution = run_candidate_analysis(
        _cli_path(inventory_path),
        _cli_path(registry_path),
        effective,
        document,
    )
    logger.info("Wrote candidate-review evidence to {}", execution.result_path)


def _registry_row(record: RegistryRecord) -> dict[str, object]:
    return {
        "database_id": record.id,
        "relative_path": record.relative_path,
        "filename": record.filename,
        "database_name": record.dbname,
        "detail_level": record.detail_level.value,
        "entries": record.entry_count,
        "targets": record.target_count,
        "decoys": record.decoy_count,
        "contaminants": record.contaminant_count,
        "entrapments": record.entrapment_count,
        "sentinels": record.sentinel_count,
        "distinct_target_ids": record.distinct_target_ids,
        "distinct_target_sequences": record.distinct_target_sequences,
        "target_id_fingerprint": record.target_id_fingerprint,
        "target_description_fingerprint": record.target_description_fingerprint,
        "target_content_fingerprint": record.target_content_fingerprint,
        "total_residues": record.total_residues,
        "indexed_at": record.indexed_at,
    }


def _require_registry_schema(connection: RegistryConnection, registry_path: Path) -> None:
    schema_version = connection.schema_version()
    if schema_version != SCHEMA_VERSION:
        raise RegistrySchemaError(schema_version, path=registry_path)


@app.command
def registry(registry_path: Path, table_path: Path) -> None:
    """Export indexed database summaries and versioned fingerprints.

    Args:
        registry_path: Existing ``.sqlite3`` or ``.duckdb`` registry.
        table_path: CSV, TSV, XLSX, or Parquet output.
    """
    with connect_registry(registry_path, read_only=True) as connection:
        _require_registry_schema(connection, registry_path)
        records = list_databases(connection)
    frame = pl.DataFrame([_registry_row(record) for record in records])
    _writer_for(table_path)(frame, table_path)
    logger.info("Wrote {} registry rows to {}", frame.height, table_path)


@app.command
def compare(
    registry_path: Path,
    database_id: int,
    table_path: Path,
    *,
    kind: EntryKind = EntryKind.TARGET,
    threshold: float = 0.99,
) -> None:
    """Compare one indexed database with every detailed database.

    Args:
        registry_path: Existing ``.sqlite3`` or ``.duckdb`` registry.
        database_id: Integer database identifier shown by ``registry``.
        table_path: CSV, TSV, XLSX, or Parquet output.
        kind: Compare target or contaminant entries.
        threshold: Containment threshold for relationship labels.
    """
    rows = [
        {
            "database_id": item.database.database_id,
            "filename": item.database.filename,
            "kind": item.kind.value,
            "relationship": item.relationship.value,
            "selected_ids": item.selected_ids,
            "other_ids": item.other_ids,
            "shared_ids": item.shared_ids,
            "id_jaccard": item.id_jaccard,
            "containment": item.containment,
            "shared_sequences": item.shared_sequence_checksums,
            "sequence_jaccard": item.sequence_jaccard,
            "shared_descriptions": item.shared_descriptions,
            "description_jaccard": item.description_jaccard,
            "shared_exact_pairs": item.shared_exact_pairs,
            "changed_shared_ids": item.changed_shared_ids,
            "exact_content": item.exact_content,
        }
        for item in compare_database(
            registry_path,
            database_id,
            threshold,
            kind=kind,
        )
    ]
    frame = pl.DataFrame(rows)
    _writer_for(table_path)(frame, table_path)
    logger.info("Wrote {} database comparisons to {}", frame.height, table_path)


@app.command
def pairs(
    registry_path: Path,
    output: Path,
    *,
    ids: Path | None = None,
    sequences: Path | None = None,
) -> None:
    """Export materialized pair metrics as stable TSV files.

    Args:
        registry_path: Existing ``.sqlite3`` or ``.duckdb`` registry.
        output: Long-form pair-statistics TSV.
        ids: Optional target-ID Jaccard matrix TSV.
        sequences: Optional target-sequence Jaccard matrix TSV.
    """
    with connect_registry(registry_path, read_only=True) as connection:
        _require_registry_schema(connection, registry_path)
        data = query_similarity_data(connection)
    write_similarity_exports(data, output, id_matrix=ids, sequence_matrix=sequences)
    logger.info("Wrote {} pair rows to {}", len(data.pairs), output)


@app.command
def cluster(
    registry_path: Path,
    table_path: Path,
    *,
    metric: ClusteringMetric = ClusteringMetric.TARGET_IDS,
) -> None:
    """Export deterministic average-linkage clustering merges.

    Args:
        registry_path: Existing ``.sqlite3`` or ``.duckdb`` registry.
        table_path: CSV, TSV, XLSX, or Parquet output.
        metric: Target identifier or target sequence Jaccard metric.
    """
    with connect_registry(registry_path, read_only=True) as connection:
        _require_registry_schema(connection, registry_path)
        result = cluster_registered_databases(connection, metric=metric)
    frame = pl.DataFrame(
        [
            {
                "metric": result.metric.value,
                "cluster_id": merge.cluster_id,
                "left_id": merge.left_id,
                "right_id": merge.right_id,
                "distance": merge.distance,
                "leaf_count": merge.leaf_count,
            }
            for merge in result.merges
        ],
        schema={
            "metric": pl.String,
            "cluster_id": pl.Int64,
            "left_id": pl.Int64,
            "right_id": pl.Int64,
            "distance": pl.Float64,
            "leaf_count": pl.Int64,
        },
    )
    _writer_for(table_path)(frame, table_path)
    logger.info("Wrote {} cluster merges to {}", frame.height, table_path)
    logger.info("Leaf order: {}", ", ".join(result.ordered_relative_paths) or "none")


def main() -> int:
    """Run the console application."""
    result = app()
    return int(result) if result is not None else 0


if __name__ == "__main__":
    sys.exit(main())
