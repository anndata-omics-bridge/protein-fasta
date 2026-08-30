"""Command-line inspection of protein FASTA databases."""

from __future__ import annotations

import datetime
import os
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Literal

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
from protein_fasta.analytics_compile import compile_digestion, make_digestion
from protein_fasta.artifact_io import write_json_atomic
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
    load_diagnostic_document,
    load_digestion_document,
    load_entry_classifier_document,
    load_enzyme_document,
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
from protein_fasta.schema.base import DocumentBase
from protein_fasta.schema.build import DatabaseBuildRequestDocument
from protein_fasta.schema.candidate import CandidateRequestDocument
from protein_fasta.schema.decoy import (
    DecoyRequestDocument,
    DecoyStrategyDocument,
    ReverseDecoyDocument,
    ShuffleDecoyDocument,
)
from protein_fasta.schema.decoy_report import DecoyReportRequestDocument
from protein_fasta.schema.peptide import (
    MemoryPeptideExecutionDocument,
    PeptideBuildRequestDocument,
    PeptideComparisonRequestDocument,
)
from protein_fasta.schema.protein_input import (
    DerivedProteinInputRequestDocument,
    ProteinInputRequestDocument,
    TargetProteinSourceDocument,
)
from protein_fasta.schema.registry import RegistryDocument
from protein_fasta.schema.uniprot import (
    CanonicalGeneDownloadDocument,
    CompleteDownloadDocument,
    ProteomeIdSelectionDocument,
    ReferenceProteomesDocument,
    ReviewedDownloadDocument,
    UniProtAcquisitionDocument,
    UniProtCatalogRequestDocument,
    UniProtDownloadRequestDocument,
)
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
type _UniProtDownloadCliMode = Literal["reviewed", "canonical", "opg"]
type _DecoyCliMethod = Literal["reverse", "shuffle"]

_UNIPROT_ACQUISITION_BY_CLI_MODE: dict[
    _UniProtDownloadCliMode,
    UniProtAcquisitionDocument,
] = {
    "reviewed": ReviewedDownloadDocument(),
    "canonical": CompleteDownloadDocument(),
    "opg": CanonicalGeneDownloadDocument(),
}
_DECOY_STRATEGY_BY_CLI_METHOD: dict[_DecoyCliMethod, DecoyStrategyDocument] = {
    "reverse": ReverseDecoyDocument(),
    "shuffle": ShuffleDecoyDocument(seed=2000),
}

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


def _recorded_path(path: Path, base: Path, /) -> Path:
    """Return a path recorded relative to its authored request when possible."""
    return Path(os.path.relpath(_cli_path(path), _cli_path(base)))


def _request_path(primary_output: Path, save: Path | None, /) -> Path:
    """Select the non-overwriting authored-request destination."""
    if save is not None:
        return _cli_path(save)
    output = _cli_path(primary_output)
    return output.with_suffix(f"{output.suffix}.request.json")


def _directory_request_path(output_dir: Path, command: str, save: Path | None, /) -> Path:
    """Select an authored-request path for a multi-artifact output directory."""
    if save is not None:
        return _cli_path(save)
    return _cli_path(output_dir) / f"{command}.request.json"


def _write_authored_request(path: Path, document: DocumentBase, /) -> None:
    """Persist one validated direct-CLI request without replacing prior intent."""
    payload = document.model_dump(mode="json", exclude_none=True)
    write_json_atomic(
        path,
        payload,
        replace_existing=False,
    )
    logger.info("Wrote authored request to {}", path)


def _reject_replay_options(request: Path | None, **options: object) -> None:
    """Reject direct-only parameters when an authored request is replayed."""
    supplied = sorted(name for name, value in options.items() if value is not None)
    if request is not None and supplied:
        rendered = ", ".join(f"--{name.replace('_', '-')}" for name in supplied)
        raise ValueError(f"{rendered} cannot be combined with --request")


def _require_direct(value: object | None, option: str, /) -> None:
    """Require one direct-mode value after replay mode has been excluded."""
    if value is None:
        raise ValueError(f"direct mode requires --{option}")


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
def diagnostics(
    fasta_path: Path,
    *,
    rules: Path | None = None,
    classifiers: Path | None = None,
) -> None:
    """Summarize built-in record diagnostics for one FASTA database.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
        rules: Optional explicit diagnostic-rules JSON.
        classifiers: Optional explicit entry-classifier JSON.
    """
    runtime_rules = make_diagnostic_rules(
        load_builtin_diagnostic_document() if rules is None else load_diagnostic_document(rules),
        load_builtin_entry_classifier_document()
        if classifiers is None
        else load_entry_classifier_document(classifiers),
    )
    summary = summarize_protein_diagnostics(iter_protein_diagnostics(fasta_path, runtime_rules))
    _log_diagnostics_summary(summary)


@app.command
def digest(
    fasta_path: Path,
    table_path: Path,
    *,
    config: Path | None = None,
    rules: Path | None = None,
) -> None:
    """Export theoretical peptides and missed-cleavage counts.

    Args:
        fasta_path: Plain, gzip, or bzip2 protein FASTA input.
        table_path: Output path ending in ``.csv``, ``.tsv``, ``.xlsx``, or ``.parquet``.
        config: Optional digestion JSON; packaged trypsin defaults are used when omitted.
        rules: Optional explicit enzyme-rules JSON.
    """
    document = load_digestion_document(config) if config is not None else DigestionDocument()
    digestion = (
        make_digestion(document)
        if rules is None
        else compile_digestion(document, load_enzyme_document(rules))
    )
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
def peptides(
    inventory_path: Path,
    *,
    request: Path | None = None,
    output: Path | None = None,
    enzyme: str | None = None,
    minimum: int | None = None,
    maximum: int | None = None,
    missed: int | None = None,
    workers: int | None = None,
    partition: int | None = None,
    save: Path | None = None,
) -> None:
    """Build peptide artifacts from direct arguments or an authored request.

    Args:
        inventory_path: Canonical biological or search inventory Parquet.
        request: Existing peptide-build request JSON to replay.
        output: Direct output directory, or replay output-directory override.
        enzyme: Direct enzyme name; packaged trypsin is the default.
        minimum: Direct minimum peptide length.
        maximum: Direct maximum peptide length.
        missed: Direct maximum missed cleavages.
        workers: Direct worker count for the memory executor.
        partition: Direct proteins per deterministic partition.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(
        request,
        enzyme=enzyme,
        minimum=minimum,
        maximum=maximum,
        missed=missed,
        workers=workers,
        partition=partition,
        save=save,
    )
    if request is None:
        _require_direct(output, "output")
        assert output is not None
        request_path = _directory_request_path(output, "peptides", save)
        request_base = request_path.parent
        destination = _cli_path(output)
        digestion = DigestionDocument(
            enzyme="trypsin" if enzyme is None else enzyme,
            min_length=7 if minimum is None else minimum,
            max_length=50 if maximum is None else maximum,
            missed_cleavages=0 if missed is None else missed,
        )
        execution = MemoryPeptideExecutionDocument(
            workers=1 if workers is None else workers,
            partition_size=500 if partition is None else partition,
        )
        document = PeptideBuildRequestDocument(
            peptides_parquet=_recorded_path(destination / "peptides.parquet", request_base),
            mapping_parquet=_recorded_path(
                destination / "protein-peptide-map.parquet",
                request_base,
            ),
            peptide_fasta=_recorded_path(destination / "peptides.fasta", request_base),
            digestion=digestion,
            execution=execution,
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_peptide_build_request(request_path)
        if output is not None:
            destination = _cli_path(output)
            document = document.model_copy(
                update={
                    "peptides_parquet": destination / document.peptides_parquet.name,
                    "mapping_parquet": destination / document.mapping_parquet.name,
                    "peptide_fasta": destination / document.peptide_fasta.name,
                }
            )
    effective = resolve_peptide_build_request(
        document,
        request_base=request_base,
    )
    execution = run_peptide_build(_cli_path(inventory_path), effective)
    logger.info("Wrote peptide-build evidence to {}", execution.result_path)


@app.command
def pepcompare(
    peptides_a: Path,
    peptides_b: Path,
    *,
    request: Path | None = None,
    output: Path | None = None,
    save: Path | None = None,
) -> None:
    """Compare peptide inventories from direct arguments or an authored request.

    Args:
        peptides_a: First canonical peptides Parquet.
        peptides_b: Second canonical peptides Parquet.
        request: Existing peptide-comparison request JSON to replay.
        output: Direct comparison Parquet, or replay output override.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(request, save=save)
    if request is None:
        _require_direct(output, "output")
        assert output is not None
        destination = _cli_path(output)
        request_path = _request_path(destination, save)
        document = PeptideComparisonRequestDocument(
            output_parquet=_recorded_path(destination, request_path.parent)
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        document = load_peptide_comparison_request(request_path)
        if output is not None:
            document = document.model_copy(update={"output_parquet": _cli_path(output)})
    effective = resolve_peptide_comparison_request(
        document,
        request_base=request_path.parent,
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
    *,
    request: Path | None = None,
    output: Path | None = None,
    timeout: float | None = None,
    save: Path | None = None,
) -> None:
    """Synchronize a UniProt catalog from direct arguments or an authored request.

    Args:
        request: Existing UniProt catalog request JSON to replay.
        output: Direct catalog directory, or replay output-directory override.
        timeout: Direct provider timeout in seconds.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(request, timeout=timeout, save=save)
    if request is None:
        _require_direct(output, "output")
        assert output is not None
        request_path = _directory_request_path(output, "uniprot-catalog", save)
        request_base = request_path.parent
        document = UniProtCatalogRequestDocument(
            output_dir=_recorded_path(output, request_base),
            selection=ReferenceProteomesDocument(),
            timeout_seconds=120.0 if timeout is None else timeout,
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_uniprot_catalog_request(request_path)
        if output is not None:
            document = document.model_copy(update={"output_dir": _cli_path(output)})
    execution = sync_uniprot_catalog(
        document,
        request_base=request_base,
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
    source: str | None = None,
    mode: _UniProtDownloadCliMode | None = None,
    *,
    request: Path | None = None,
    output: Path | None = None,
    timeout: float | None = None,
    save: Path | None = None,
) -> None:
    """Acquire one UniProt proteome FASTA from direct arguments or a request document.

    Args:
        source: Direct UniProt proteome identifier.
        mode: Direct mode: reviewed, canonical, or one protein per gene (opg).
        request: Existing UniProt download request JSON to replay.
        output: Direct FASTA destination, or replay output override.
        timeout: Direct provider timeout in seconds.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(request, source=source, mode=mode, timeout=timeout, save=save)
    if request is None:
        _require_direct(source, "source")
        _require_direct(mode, "mode")
        assert source is not None
        assert mode is not None
        destination = _cli_path(output if output is not None else Path(f"{source}_{mode}.fasta"))
        request_path = _request_path(destination, save)
        document = UniProtDownloadRequestDocument(
            selection=ProteomeIdSelectionDocument(proteome_id=source),
            acquisition=_UNIPROT_ACQUISITION_BY_CLI_MODE[mode],
            output_fasta=_recorded_path(destination, request_path.parent),
            timeout_seconds=120.0 if timeout is None else timeout,
        )
        _write_authored_request(request_path, document)
        request_base = request_path.parent
        output_override = None
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_uniprot_download_request(request_path)
        output_override = None if output is None else _cli_path(output)
    effective = resolve_uniprot_download(
        document,
        request_base=request_base,
        overrides=UniProtDownloadOverrides(
            output_fasta=output_override,
            timeout_seconds=None,
        ),
    )
    execution = run_uniprot_download(effective)
    logger.info("Wrote UniProt download evidence to {}", execution.result_path)


@app.command
def prepare(
    fasta_path: Path | None = None,
    output: Path | None = None,
    *,
    id: str | None = None,
    request: Path | None = None,
    save: Path | None = None,
) -> None:
    """Prepare one target FASTA directly or replay an authored source request.

    Args:
        fasta_path: Direct target FASTA source.
        output: Direct protein-input Parquet, or replay output override.
        id: Direct stable source identifier.
        request: Existing protein-input request JSON to replay.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(request, fasta_path=fasta_path, id=id, save=save)
    if request is None:
        _require_direct(fasta_path, "fasta-path")
        _require_direct(output, "output")
        _require_direct(id, "id")
        assert fasta_path is not None
        assert output is not None
        assert id is not None
        destination = _cli_path(output)
        request_path = _request_path(destination, save)
        request_base = request_path.parent
        document = ProteinInputRequestDocument(
            sources=(
                TargetProteinSourceDocument(
                    source_id=id,
                    path=_recorded_path(fasta_path, request_base),
                ),
            ),
            output_parquet=_recorded_path(destination, request_base),
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_protein_input_request(request_path)
        if output is not None:
            document = document.model_copy(update={"output_parquet": _cli_path(output)})
    effective = resolve_protein_input_request(
        document,
        request_base=request_base,
    )
    execution = run_protein_input_preparation(effective)
    logger.info("Wrote protein-input evidence to {}", execution.result_path)


@app.command(name="derive-input")
def derive_input(
    source_inventory: Path | None = None,
    output: Path | None = None,
    *,
    id: str | None = None,
    request: Path | None = None,
    save: Path | None = None,
) -> None:
    """Derive clean source rows directly or replay an authored request.

    The source retains target and contaminant proteins. Existing sentinel,
    section-marker, entrapment, and decoy rows are excluded. An optional second
    inventory supplies foreign proteins for a subsequent entrapment build.

    Args:
        source_inventory: Direct biological or search inventory source.
        output: Direct protein-input Parquet, or replay output override.
        id: Direct stable source identifier.
        request: Existing derived protein-input request JSON to replay.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(request, source_inventory=source_inventory, id=id, save=save)
    if request is None:
        _require_direct(source_inventory, "source-inventory")
        _require_direct(output, "output")
        _require_direct(id, "id")
        assert source_inventory is not None
        assert output is not None
        assert id is not None
        destination = _cli_path(output)
        request_path = _request_path(destination, save)
        request_base = request_path.parent
        document = DerivedProteinInputRequestDocument(
            source_inventory=_recorded_path(source_inventory, request_base),
            source_id=id,
            output_parquet=_recorded_path(destination, request_base),
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_derived_protein_input_request(request_path)
        if output is not None:
            document = document.model_copy(update={"output_parquet": _cli_path(output)})
    effective = resolve_derived_protein_input_request(
        document,
        request_base=request_base,
    )
    execution = run_derived_protein_input_preparation(effective)
    logger.info("Wrote derived protein-input evidence to {}", execution.result_path)


@app.command
def decoy(
    biological_inventory: Path,
    *,
    request: Path | None = None,
    output: Path | None = None,
    method: _DecoyCliMethod | None = None,
    prefix: str | None = None,
    save: Path | None = None,
) -> None:
    """Generate a search database from direct arguments or an authored request.

    Args:
        biological_inventory: Existing decoy-free protein-inventory Parquet.
        request: Existing decoy request JSON to replay.
        output: Direct search FASTA, or replay output override.
        method: Direct decoy method.
        prefix: Direct generated-header prefix.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(request, method=method, prefix=prefix, save=save)
    if request is None:
        _require_direct(output, "output")
        _require_direct(method, "method")
        assert output is not None
        assert method is not None
        destination = _cli_path(output)
        request_path = _request_path(destination, save)
        request_base = request_path.parent
        document = DecoyRequestDocument(
            output_fasta=_recorded_path(destination, request_base),
            decoy_prefix="REV_" if prefix is None else prefix,
            strategy=_DECOY_STRATEGY_BY_CLI_METHOD[method],
        )
        _write_authored_request(request_path, document)
        output_override = None
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_decoy_request(request_path)
        output_override = None if output is None else _cli_path(output)
    effective = resolve_decoy_request(
        document,
        request_base=request_base,
        overrides=DecoyOverrides(
            output_fasta=output_override,
            decoy_prefix=None,
        ),
    )
    execution = run_decoy_generation(_cli_path(biological_inventory), effective)
    logger.info("Wrote decoy-generation evidence to {}", execution.result_path)


@app.command(name="decoy-report")
def decoy_report(
    biological_inventory: Path,
    *,
    request: Path | None = None,
    output: Path | None = None,
    method: tuple[_DecoyCliMethod, ...] = (),
    prefix: str | None = None,
    enzyme: str | None = None,
    minimum: int | None = None,
    maximum: int | None = None,
    missed: int | None = None,
    save: Path | None = None,
) -> None:
    """Compare decoy methods from direct arguments or an authored request.

    Args:
        biological_inventory: Existing decoy-free protein-inventory Parquet.
        request: Existing decoy-method comparison request JSON to replay.
        output: Direct comparison Parquet, or replay output override.
        method: One or more direct decoy methods.
        prefix: Direct generated-header prefix.
        enzyme: Direct collision-digestion enzyme.
        minimum: Direct minimum peptide length.
        maximum: Direct maximum peptide length.
        missed: Direct maximum missed cleavages.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(
        request,
        method=method or None,
        prefix=prefix,
        enzyme=enzyme,
        minimum=minimum,
        maximum=maximum,
        missed=missed,
        save=save,
    )
    if request is None:
        _require_direct(output, "output")
        if not method:
            raise ValueError("direct mode requires at least one --method")
        assert output is not None
        destination = _cli_path(output)
        request_path = _request_path(destination, save)
        request_base = request_path.parent
        document = DecoyReportRequestDocument(
            output_parquet=_recorded_path(destination, request_base),
            decoy_prefix="REV_" if prefix is None else prefix,
            digestion=DigestionDocument(
                enzyme="trypsin" if enzyme is None else enzyme,
                min_length=7 if minimum is None else minimum,
                max_length=50 if maximum is None else maximum,
                missed_cleavages=0 if missed is None else missed,
            ),
            strategies=tuple(_DECOY_STRATEGY_BY_CLI_METHOD[name] for name in method),
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_decoy_report_request(request_path)
        if output is not None:
            document = document.model_copy(update={"output_parquet": _cli_path(output)})
    effective = resolve_decoy_report_request(
        document,
        request_base=request_base,
    )
    execution = run_decoy_report(_cli_path(biological_inventory), effective)
    logger.info("Wrote decoy-method evidence to {}", execution.result_path)


@app.command
def build(
    protein_input: Path,
    *,
    request: Path | None = None,
    output: Path | None = None,
    project: int | None = None,
    dbn: int | None = None,
    description: str | None = None,
    profile: Path | None = None,
    date: datetime.date | None = None,
    save: Path | None = None,
) -> None:
    """Build a biological FASTA from direct arguments or an authored request.

    Relative request paths resolve beside the request file. Profile defaults are
    overridden by request values and then by explicitly supplied CLI options. The
    effective request is written before sequence work starts; the final result JSON
    records artifacts, checksums, counts, summaries, normalization, and entrapment.
    Decoy generation is the separate ``decoy`` command.

    Args:
        protein_input: Canonical decoy-free protein-input Parquet.
        request: Existing biological-build request JSON to replay.
        output: Direct build directory, or replay output-directory override.
        project: Direct project identifier for packaged FGCZ naming.
        dbn: Direct project database number for packaged FGCZ naming.
        description: Direct compact database-description token.
        profile: Optional portable defaults JSON; packaged FGCZ defaults are used otherwise.
        date: Direct build date; today is used when omitted.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(
        request,
        project=project,
        dbn=dbn,
        description=description,
        date=date,
        save=save,
    )
    if request is None:
        _require_direct(output, "output")
        _require_direct(project, "project")
        _require_direct(dbn, "dbn")
        _require_direct(description, "description")
        assert output is not None
        assert project is not None
        assert dbn is not None
        assert description is not None
        request_path = _directory_request_path(output, "build", save)
        request_base = request_path.parent
        document = DatabaseBuildRequestDocument(
            output_dir=_recorded_path(output, request_base),
            date=datetime.date.today() if date is None else date,
            name_fields={"project": project, "dbn": dbn, "description": description},
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_database_build_request(request_path)
        if output is not None:
            document = document.model_copy(update={"output_dir": _cli_path(output)})
    if profile is None:
        build_profile = load_builtin_database_build_profile()
        profile_base = request_base
    else:
        build_profile = load_database_build_profile(profile)
        profile_base = _cli_path(profile).parent
    effective = resolve_database_build(
        build_profile,
        document,
        profile_base=profile_base,
        request_base=request_base,
        overrides=DatabaseBuildOverrides(date=None),
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
    *,
    request: Path | None = None,
    output: Path | None = None,
    threshold: float | None = None,
    metric: Literal["target_ids", "target_sequences"] | None = None,
    limit: int | None = None,
    config: Path | None = None,
    save: Path | None = None,
) -> None:
    """Review a candidate from direct arguments or an authored request.

    Args:
        inventory_path: Canonical protein- or search-inventory Parquet.
        registry_path: Existing SQLite or DuckDB registry.
        request: Existing candidate-review request JSON to replay.
        output: Direct comparison Parquet, or replay output override.
        threshold: Direct relationship-overlap threshold.
        metric: Direct clustering metric.
        limit: Direct nearest-neighbour limit.
        config: Optional registry-policy JSON matching the indexed registry.
        save: Optional authored-request destination for a direct run.
    """
    _reject_replay_options(
        request,
        threshold=threshold,
        metric=metric,
        limit=limit,
        save=save,
    )
    if request is None:
        _require_direct(output, "output")
        assert output is not None
        destination = _cli_path(output)
        request_path = _request_path(destination, save)
        request_base = request_path.parent
        document = CandidateRequestDocument(
            output_parquet=_recorded_path(destination, request_base),
            overlap_threshold=0.99 if threshold is None else threshold,
            clustering_metric="target_ids" if metric is None else metric,
            neighbour_limit=50 if limit is None else limit,
        )
        _write_authored_request(request_path, document)
    else:
        request_path = _cli_path(request)
        request_base = request_path.parent
        document = load_candidate_request(request_path)
        if output is not None:
            document = document.model_copy(update={"output_parquet": _cli_path(output)})
    effective = resolve_candidate_request(
        document,
        request_base=request_base,
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
