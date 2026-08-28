"""Compose one reproducible UniProt proteome FASTA acquisition."""

from __future__ import annotations

import importlib.metadata
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from protein_fasta.artifact_io import (
    artifact_document,
    publish_exclusive,
    temporary_sibling,
    write_json_atomic,
)
from protein_fasta.schema.uniprot import (
    CompleteDownloadDocument,
    EffectiveUniProtDownloadDocument,
    ResolvedProteomeDocument,
    ReviewedDownloadDocument,
    TaxonomySelectionDocument,
    UniProtDownloadRequestDocument,
    UniProtDownloadResultDocument,
)
from protein_fasta.uniprot.acquisition import (
    CanonicalGeneProteins,
    CompleteProteins,
    ProteinAcquisition,
    ReviewedProteins,
)
from protein_fasta.uniprot.models import UniProtDownloadOutcome
from protein_fasta.uniprot.resolution import (
    ProteomeResolution,
    ResolveProteomeId,
    ResolveTaxonomy,
)
from protein_fasta.uniprot.transport import UniProtTransport


@dataclass(frozen=True, slots=True)
class UniProtDownloadOverrides:
    """Common explicitly supplied values taking precedence over one request."""

    output_fasta: Path | None = None
    timeout_seconds: float | None = None


DEFAULT_UNIPROT_DOWNLOAD_OVERRIDES = UniProtDownloadOverrides()


@dataclass(frozen=True, slots=True)
class UniProtDownloadPlan:
    """Compiled resolution and acquisition behaviors for one provider request."""

    resolution: ProteomeResolution
    acquisition: ProteinAcquisition


@dataclass(frozen=True, slots=True)
class UniProtDownloadExecution:
    """Runtime outcome and durable evidence for one completed acquisition."""

    outcome: UniProtDownloadOutcome
    document: UniProtDownloadResultDocument
    fasta_path: Path
    effective_request_path: Path
    result_path: Path


def resolve_uniprot_download(
    request: UniProtDownloadRequestDocument,
    /,
    *,
    request_base: Path,
    overrides: UniProtDownloadOverrides = DEFAULT_UNIPROT_DOWNLOAD_OVERRIDES,
) -> EffectiveUniProtDownloadDocument:
    """Resolve request-relative paths and explicit common overrides."""
    output = overrides.output_fasta or request.output_fasta
    if not output.is_absolute():
        output = request_base / output
    return EffectiveUniProtDownloadDocument(
        selection=request.selection,
        acquisition=request.acquisition,
        output_fasta=output.resolve(),
        timeout_seconds=(
            overrides.timeout_seconds
            if overrides.timeout_seconds is not None
            else request.timeout_seconds
        ),
    )


def make_uniprot_download_plan(
    effective: EffectiveUniProtDownloadDocument,
    /,
) -> UniProtDownloadPlan:
    """Compile passive selection and acquisition variants exactly once."""
    selection = effective.selection
    if isinstance(selection, TaxonomySelectionDocument):
        resolution: ProteomeResolution = ResolveTaxonomy(selection.taxid)
    else:
        resolution = ResolveProteomeId(selection.proteome_id)

    acquisition_document = effective.acquisition
    if isinstance(acquisition_document, ReviewedDownloadDocument):
        acquisition: ProteinAcquisition = ReviewedProteins()
    elif isinstance(acquisition_document, CompleteDownloadDocument):
        acquisition = CompleteProteins()
    else:
        acquisition = CanonicalGeneProteins()
    return UniProtDownloadPlan(resolution, acquisition)


def run_uniprot_download(
    effective: EffectiveUniProtDownloadDocument,
    /,
    *,
    transport: UniProtTransport | None = None,
) -> UniProtDownloadExecution:
    """Acquire, validate, and atomically publish one UniProt FASTA and evidence."""
    output = effective.output_fasta
    result_path = output.with_suffix(f"{output.suffix}.result.json")
    effective_path = output.with_suffix(f"{output.suffix}.effective.json")
    _refuse_existing(output, result_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(
        effective_path,
        effective.model_dump(mode="json"),
        replace_existing=True,
    )
    plan = make_uniprot_download_plan(effective)
    if transport is None:
        with UniProtTransport(timeout_seconds=effective.timeout_seconds) as owned_transport:
            return _run_download(plan, effective, effective_path, result_path, owned_transport)
    return _run_download(plan, effective, effective_path, result_path, transport)


def _run_download(
    plan: UniProtDownloadPlan,
    effective: EffectiveUniProtDownloadDocument,
    effective_path: Path,
    result_path: Path,
    transport: UniProtTransport,
) -> UniProtDownloadExecution:
    output = effective.output_fasta
    with temporary_sibling(output) as staged:
        proteome = plan.resolution.resolve(transport)
        outcome = plan.acquisition.acquire(transport, proteome, staged)
        _validate_download(outcome, staged)
        artifact = artifact_document(
            staged,
            recorded_path=Path(output.name),
            schema_name="uniprot-fasta",
            schema_version="1",
            row_count=outcome.transfer.actual_entry_count,
        )
        document = UniProtDownloadResultDocument(
            protein_fasta_version=_package_version(),
            effective_request=effective,
            resolved_proteome=ResolvedProteomeDocument(
                proteome_id=proteome.proteome_id,
                taxid=proteome.taxid,
                protein_count=proteome.protein_count,
                gene_count=proteome.gene_count,
                organism=proteome.organism,
                resolution_method=proteome.resolution_method.value,
                resolution_query=proteome.resolution_query,
            ),
            provider_query=outcome.provider_query,
            observed_releases=outcome.transfer.observed_releases,
            actual_entry_count=outcome.transfer.actual_entry_count,
            provider_reported_counts=outcome.transfer.provider_reported_counts,
            artifact=artifact,
            warnings=_transfer_warnings(outcome),
        )
        publish_exclusive(staged, output)
        try:
            write_json_atomic(
                result_path,
                document.model_dump(mode="json"),
                replace_existing=False,
            )
        except BaseException:
            output.unlink(missing_ok=True)
            raise
    logger.info(
        "Downloaded {} UniProt records ({}) to {}",
        outcome.transfer.actual_entry_count,
        outcome.mode.value,
        output,
    )
    return UniProtDownloadExecution(outcome, document, output, effective_path, result_path)


def _validate_download(outcome: UniProtDownloadOutcome, staged: Path) -> None:
    if outcome.transfer.actual_entry_count <= 0:
        raise ValueError(
            f"UniProt returned no FASTA entries for {outcome.resolved_proteome.proteome_id} "
            f"({outcome.mode.value})"
        )
    if not staged.is_file() or staged.stat().st_size == 0:
        raise ValueError(
            f"UniProt returned an empty FASTA for {outcome.resolved_proteome.proteome_id} "
            f"({outcome.mode.value})"
        )


def _transfer_warnings(outcome: UniProtDownloadOutcome) -> tuple[str, ...]:
    transfer = outcome.transfer
    warnings: list[str] = []
    if len(transfer.observed_releases) > 1:
        warnings.append(
            "UniProt pages reported different releases: " + ", ".join(transfer.observed_releases)
        )
    if len(transfer.provider_reported_counts) > 1:
        warnings.append(
            "UniProt responses reported different total counts: "
            + ", ".join(str(value) for value in transfer.provider_reported_counts)
        )
    if transfer.provider_reported_counts and any(
        count != transfer.actual_entry_count for count in transfer.provider_reported_counts
    ):
        warnings.append(
            f"UniProt reported {transfer.provider_reported_counts!r} records but "
            f"{transfer.actual_entry_count} FASTA records were written"
        )
    return tuple(warnings)


def _refuse_existing(output: Path, result_path: Path) -> None:
    existing = [path for path in (output, result_path) if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace existing UniProt artifacts: {names}")


def _package_version() -> str:
    return importlib.metadata.version("protein-fasta")
