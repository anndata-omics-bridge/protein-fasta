"""Synchronize, read, and filter reproducible UniProt proteome catalogs."""

from __future__ import annotations

import datetime
import importlib.metadata
import json
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import cast

import polars as pl
from loguru import logger
from pydantic import ValidationError

from protein_fasta.analytics.hashing import file_checksum
from protein_fasta.artifact_io import (
    artifact_document,
    publish_exclusive,
    temporary_sibling,
    write_json_atomic,
)
from protein_fasta.schema.uniprot import (
    AllProteomesDocument,
    ReferenceProteomesDocument,
    UniProtCatalogRequestDocument,
    UniProtCatalogResultDocument,
)
from protein_fasta.uniprot.models import UniProtCatalogRow
from protein_fasta.uniprot.provider_rows import catalog_row
from protein_fasta.uniprot.transport import UniProtTransport

_CATALOG_PREFIX = "uniprot-proteomes-"
_CATALOG_RESULT_PREFIX = "uniprot-catalog-"
_TIMESTAMP_FORMAT = "%Y%m%dT%H%M%S%fZ"

_CATALOG_SCHEMA = pl.Schema(
    {
        "proteome_id": pl.String,
        "taxid": pl.Int64,
        "organism": pl.String,
        "proteome_type": pl.String,
        "swissprot": pl.Int64,
        "swissprot_trembl": pl.Int64,
        "one_seq_per_gene": pl.Int64,
    }
)

type CatalogProgress = Callable[[int], None]


@dataclass(frozen=True, slots=True)
class UniProtCatalogExecution:
    """Canonical frame and durable artifacts from one catalog synchronization."""

    frame: pl.DataFrame
    document: UniProtCatalogResultDocument
    catalog_path: Path
    result_path: Path


def ignore_catalog_progress(count: int, /) -> None:
    """Identity progress sink for callers that do not render progress."""


def catalog_query(request: UniProtCatalogRequestDocument, /) -> str:
    """Compile one catalog selection variant into its exact provider query."""
    selection = request.selection
    if isinstance(selection, ReferenceProteomesDocument):
        return "reference:true"
    if isinstance(selection, AllProteomesDocument):
        return "*"
    return selection.query


def sync_uniprot_catalog(
    request: UniProtCatalogRequestDocument,
    /,
    *,
    request_base: Path,
    transport: UniProtTransport | None = None,
    retrieved_at: datetime.datetime | None = None,
    progress: CatalogProgress = ignore_catalog_progress,
) -> UniProtCatalogExecution:
    """Fetch and atomically publish one immutable UniProt proteome catalog."""
    moment = retrieved_at or datetime.datetime.now(datetime.UTC)
    if moment.tzinfo is None:
        raise ValueError("retrieved_at must be timezone-aware")
    moment = moment.astimezone(datetime.UTC)
    stamp = moment.strftime(_TIMESTAMP_FORMAT)
    output_dir = request.output_dir
    if not output_dir.is_absolute():
        output_dir = request_base / output_dir
    output_dir = output_dir.resolve()
    effective_request = request.model_copy(update={"output_dir": output_dir})
    catalog_path = output_dir / f"{_CATALOG_PREFIX}{stamp}.parquet"
    result_path = output_dir / f"{_CATALOG_RESULT_PREFIX}{stamp}.result.json"
    _refuse_catalog_collision(catalog_path, result_path)
    query = catalog_query(request)
    if transport is None:
        with UniProtTransport(timeout_seconds=request.timeout_seconds) as owned_transport:
            return _sync_catalog(
                effective_request,
                query,
                moment,
                catalog_path,
                result_path,
                owned_transport,
                progress,
            )
    return _sync_catalog(
        effective_request,
        query,
        moment,
        catalog_path,
        result_path,
        transport,
        progress,
    )


def read_uniprot_catalog(path: Path, /) -> pl.DataFrame:
    """Read one canonical local catalog without contacting UniProt."""
    frame = pl.read_parquet(path)
    if frame.schema != _CATALOG_SCHEMA:
        raise ValueError(f"invalid UniProt catalog schema in {path}: {frame.schema!r}")
    return frame


def filter_uniprot_catalog(frame: pl.DataFrame, query: str, /) -> pl.DataFrame:
    """Filter a canonical catalog by proteome id, organism, or taxid text."""
    _validate_catalog_frame(frame)
    needle = query.strip().lower()
    if not needle:
        return frame
    proteome_ids = cast("list[str]", frame.get_column("proteome_id").to_list())
    organisms = cast("list[str | None]", frame.get_column("organism").to_list())
    taxids = cast("list[int | None]", frame.get_column("taxid").to_list())
    indices = [
        index
        for index, (proteome_id, organism, taxid) in enumerate(
            zip(proteome_ids, organisms, taxids, strict=True)
        )
        if needle in f"{proteome_id} {organism or ''} {taxid or ''}".lower()
    ]
    return frame[indices]


def latest_uniprot_catalog(directory: Path, /) -> Path | None:
    """Return the newest manifest-committed catalog with matching exact bytes."""
    catalogs = list_uniprot_catalogs(directory)
    return catalogs[-1] if catalogs else None


def list_uniprot_catalogs(directory: Path, /) -> list[Path]:
    """Return manifest-committed catalogs in chronological filename order."""
    if not directory.is_dir():
        return []
    catalogs: list[Path] = []
    for result_path in sorted(directory.glob(f"{_CATALOG_RESULT_PREFIX}*.result.json")):
        document = _read_catalog_result(result_path)
        catalog_path = result_path.parent / document.artifact.path
        if catalog_path.is_file() and file_checksum(catalog_path) == document.artifact.checksum:
            catalogs.append(catalog_path)
    return catalogs


def _sync_catalog(
    request: UniProtCatalogRequestDocument,
    query: str,
    moment: datetime.datetime,
    catalog_path: Path,
    result_path: Path,
    transport: UniProtTransport,
    progress: CatalogProgress,
) -> UniProtCatalogExecution:
    rows: list[UniProtCatalogRow] = []
    releases: list[str] = []
    reported_counts: list[int] = []
    for page in transport.iter_proteome_pages(query):
        if page.release is not None:
            releases.append(page.release)
        if page.reported_count is not None:
            reported_counts.append(page.reported_count)
        for raw in page.records:
            rows.append(catalog_row(raw))
            if len(rows) % 2000 == 0:
                progress(len(rows))
                logger.info("Retrieved {} UniProt proteomes", len(rows))

    frame = _catalog_frame(rows)
    with temporary_sibling(catalog_path) as staged_catalog:
        frame.write_parquet(staged_catalog)
        artifact = artifact_document(
            staged_catalog,
            recorded_path=Path(catalog_path.name),
            schema_name="uniprot-proteome-catalog",
            schema_version="0.1",
            row_count=frame.height,
        )
        document = UniProtCatalogResultDocument(
            protein_fasta_version=importlib.metadata.version("protein-fasta"),
            request=request,
            provider_query=query,
            retrieved_at=moment,
            observed_releases=tuple(dict.fromkeys(releases)),
            provider_reported_counts=tuple(dict.fromkeys(reported_counts)),
            artifact=artifact,
            warnings=_catalog_warnings(frame.height, releases, reported_counts),
        )
        publish_exclusive(staged_catalog, catalog_path)
        try:
            write_json_atomic(
                result_path,
                document.model_dump(mode="json"),
                replace_existing=False,
            )
        except BaseException:
            catalog_path.unlink(missing_ok=True)
            raise
    logger.info("Synchronized {} UniProt proteomes to {}", frame.height, catalog_path)
    return UniProtCatalogExecution(frame, document, catalog_path, result_path)


def _catalog_frame(rows: list[UniProtCatalogRow]) -> pl.DataFrame:
    return pl.DataFrame((asdict(row) for row in rows), schema=_CATALOG_SCHEMA)


def _validate_catalog_frame(frame: pl.DataFrame) -> None:
    if frame.schema != _CATALOG_SCHEMA:
        raise ValueError(f"invalid UniProt catalog frame schema: {frame.schema!r}")


def _catalog_warnings(
    actual_count: int,
    releases: list[str],
    reported_counts: list[int],
) -> tuple[str, ...]:
    distinct_releases = tuple(dict.fromkeys(releases))
    distinct_counts = tuple(dict.fromkeys(reported_counts))
    warnings: list[str] = []
    if len(distinct_releases) > 1:
        warnings.append(
            "UniProt pages reported different releases: " + ", ".join(distinct_releases)
        )
    if len(distinct_counts) > 1:
        warnings.append(
            "UniProt pages reported different total counts: "
            + ", ".join(str(value) for value in distinct_counts)
        )
    if distinct_counts and any(count != actual_count for count in distinct_counts):
        warnings.append(
            f"UniProt reported {distinct_counts!r} proteomes but {actual_count} rows were written"
        )
    return tuple(warnings)


def _read_catalog_result(path: Path) -> UniProtCatalogResultDocument:
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return UniProtCatalogResultDocument.model_validate(value)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValidationError) as error:
        raise ValueError(f"cannot read UniProt catalog result {path}: {error}") from error


def _refuse_catalog_collision(catalog_path: Path, result_path: Path) -> None:
    existing = [path for path in (catalog_path, result_path) if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to replace immutable UniProt catalog artifacts: {names}")
