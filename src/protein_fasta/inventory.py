"""Canonical biological and search inventory frame contracts."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import polars as pl

from protein_fasta.database.models import (
    BiologicalDatabase,
    BiologicalKind,
    DecoyInventoryEntry,
    DecoySourceKind,
    ProteinInventoryEntry,
    SearchDatabase,
    SearchInventoryEntry,
    SourceRole,
)

PROTEIN_INPUT_SCHEMA = pl.Schema(
    {
        "source_order": pl.Int64,
        "record_order": pl.Int64,
        "source_id": pl.String,
        "role": pl.String,
        "block_name": pl.String,
        "block_description": pl.String,
        "raw_header": pl.String,
        "id": pl.String,
        "description": pl.String,
        "sequence": pl.String,
        "upper_cased": pl.Boolean,
        "terminal_stop_stripped": pl.Boolean,
    }
)

PROTEIN_INVENTORY_SCHEMA = pl.Schema(
    {
        "final_order": pl.Int64,
        "source_order": pl.Int64,
        "record_order": pl.Int64,
        "source_id": pl.String,
        "source_role": pl.String,
        "raw_header": pl.String,
        "id": pl.String,
        "description": pl.String,
        "sequence": pl.String,
        "kind": pl.String,
        "contaminant_group": pl.String,
        "sequence_hash": pl.String,
        "entrapment_strategy": pl.String,
    }
)

SEARCH_INVENTORY_SCHEMA = pl.Schema(
    {
        **dict(PROTEIN_INVENTORY_SCHEMA),
        "decoy_strategy": pl.String,
        "decoy_source_order": pl.Int64,
        "decoy_source_id": pl.String,
        "decoy_source_kind": pl.String,
    }
)

ENTRAPMENT_PAIR_SCHEMA = pl.Schema(
    {
        "source_id": pl.String,
        "target_peptide": pl.String,
        "generated_peptide": pl.String,
        "fold_index": pl.Int64,
    }
)


def protein_inventory_frame(rows: list[dict[str, object]], /) -> pl.DataFrame:
    """Construct one canonical biological inventory frame."""
    return pl.DataFrame(rows, schema=PROTEIN_INVENTORY_SCHEMA)


def biological_database_frame(database: BiologicalDatabase, /) -> pl.DataFrame:
    """Project one runtime biological database to its canonical frame."""
    return protein_inventory_frame([_protein_row(entry) for entry in database.entries])


def read_protein_input(path: Path, /) -> pl.DataFrame:
    """Read a canonical protein input and reject schema drift."""
    frame = pl.read_parquet(path)
    validate_protein_input_frame(frame, source=str(path))
    return frame


def validate_protein_input_frame(
    frame: pl.DataFrame,
    /,
    *,
    source: str = "in-memory frame",
) -> None:
    """Reject schema drift at the in-memory protein-input boundary."""
    if frame.schema != PROTEIN_INPUT_SCHEMA:
        raise ValueError(f"invalid protein-input schema in {source}: {frame.schema!r}")


def search_inventory_frame(rows: list[dict[str, object]], /) -> pl.DataFrame:
    """Construct one canonical search inventory frame."""
    return pl.DataFrame(rows, schema=SEARCH_INVENTORY_SCHEMA)


def search_database_frame(database: SearchDatabase, /) -> pl.DataFrame:
    """Project one runtime search database to its canonical frame."""
    rows = [
        _decoy_row(entry) if isinstance(entry, DecoyInventoryEntry) else _search_protein_row(entry)
        for entry in database.entries
    ]
    return search_inventory_frame(rows)


def entrapment_pair_frame(rows: list[dict[str, object]], /) -> pl.DataFrame:
    """Construct one canonical entrapment-pair frame."""
    return pl.DataFrame(rows, schema=ENTRAPMENT_PAIR_SCHEMA)


def read_entrapment_pairs(path: Path, /) -> pl.DataFrame:
    """Read and validate one canonical entrapment-pair artifact."""
    frame = pl.read_parquet(path)
    if frame.schema != ENTRAPMENT_PAIR_SCHEMA:
        raise ValueError(f"invalid entrapment-pair schema in {path}: {frame.schema!r}")
    return frame


def read_protein_inventory(path: Path, /) -> pl.DataFrame:
    """Read a canonical biological inventory and reject schema drift."""
    frame = pl.read_parquet(path)
    if frame.schema != PROTEIN_INVENTORY_SCHEMA:
        raise ValueError(f"invalid protein-inventory schema in {path}: {frame.schema!r}")
    return frame


def read_search_inventory(path: Path, /) -> pl.DataFrame:
    """Read a canonical search inventory and reject schema drift."""
    frame = pl.read_parquet(path)
    if frame.schema != SEARCH_INVENTORY_SCHEMA:
        raise ValueError(f"invalid search-inventory schema in {path}: {frame.schema!r}")
    return frame


def read_database_inventory(path: Path, /) -> tuple[SearchInventoryEntry, ...]:
    """Read either canonical biological or search inventory into runtime entries."""
    frame = pl.read_parquet(path)
    if frame.schema == PROTEIN_INVENTORY_SCHEMA:
        return protein_inventory_entries(frame)
    if frame.schema == SEARCH_INVENTORY_SCHEMA:
        return search_inventory_entries(frame)
    raise ValueError(f"invalid biological/search inventory schema in {path}: {frame.schema!r}")


def protein_inventory_entries(frame: pl.DataFrame, /) -> tuple[ProteinInventoryEntry, ...]:
    """Project a validated biological frame once into frozen runtime entries."""
    _validate_frame(frame, PROTEIN_INVENTORY_SCHEMA, "protein inventory")
    return tuple(_protein_entry(row) for row in frame.iter_rows(named=True))


def search_inventory_entries(frame: pl.DataFrame, /) -> tuple[SearchInventoryEntry, ...]:
    """Project a validated search frame once into its runtime sum type."""
    _validate_frame(frame, SEARCH_INVENTORY_SCHEMA, "search inventory")
    entries: list[SearchInventoryEntry] = []
    for row in frame.iter_rows(named=True):
        if row["kind"] == "decoy":
            entries.append(_decoy_entry(row))
        else:
            entries.append(_protein_entry(row))
    return tuple(entries)


def _protein_entry(row: dict[str, object]) -> ProteinInventoryEntry:
    kind = str(row["kind"])
    if kind not in {"sentinel", "target", "contaminant", "entrapment"}:
        raise ValueError(f"invalid biological inventory kind: {kind!r}")
    source_role = _optional_string(row["source_role"])
    if source_role not in {None, "target", "contaminant", "foreign", "entrapment"}:
        raise ValueError(f"invalid protein source role: {source_role!r}")
    return ProteinInventoryEntry(
        final_order=int(str(row["final_order"])),
        raw_header=str(row["raw_header"]),
        identifier=str(row["id"]),
        description=_optional_string(row["description"]),
        sequence=str(row["sequence"]),
        kind=cast("BiologicalKind", kind),
        contaminant_group=_optional_string(row["contaminant_group"]),
        sequence_hash=str(row["sequence_hash"]),
        entrapment_strategy=_optional_string(row["entrapment_strategy"]),
        source_order=_optional_int(row["source_order"]),
        record_order=_optional_int(row["record_order"]),
        source_id=_optional_string(row["source_id"]),
        source_role=cast("SourceRole | None", source_role),
    )


def _decoy_entry(row: dict[str, object]) -> DecoyInventoryEntry:
    source_kind = str(row["decoy_source_kind"])
    if source_kind not in {"target", "contaminant", "entrapment"}:
        raise ValueError(f"invalid decoy source kind: {source_kind!r}")
    return DecoyInventoryEntry(
        final_order=int(str(row["final_order"])),
        raw_header=str(row["raw_header"]),
        identifier=str(row["id"]),
        description=_optional_string(row["description"]),
        sequence=str(row["sequence"]),
        contaminant_group=_optional_string(row["contaminant_group"]),
        sequence_hash=str(row["sequence_hash"]),
        entrapment_strategy=_optional_string(row["entrapment_strategy"]),
        decoy_strategy=str(row["decoy_strategy"]),
        decoy_source_order=int(str(row["decoy_source_order"])),
        decoy_source_id=str(row["decoy_source_id"]),
        decoy_source_kind=cast("DecoySourceKind", source_kind),
    )


def _protein_row(entry: ProteinInventoryEntry) -> dict[str, object]:
    return {
        "final_order": entry.final_order,
        "source_order": entry.source_order,
        "record_order": entry.record_order,
        "source_id": entry.source_id,
        "source_role": entry.source_role,
        "raw_header": entry.raw_header,
        "id": entry.identifier,
        "description": entry.description,
        "sequence": entry.sequence,
        "kind": entry.kind,
        "contaminant_group": entry.contaminant_group,
        "sequence_hash": entry.sequence_hash,
        "entrapment_strategy": entry.entrapment_strategy,
    }


def _search_protein_row(entry: ProteinInventoryEntry) -> dict[str, object]:
    return {
        **_protein_row(entry),
        "decoy_strategy": None,
        "decoy_source_order": None,
        "decoy_source_id": None,
        "decoy_source_kind": None,
    }


def _decoy_row(entry: DecoyInventoryEntry) -> dict[str, object]:
    return {
        "final_order": entry.final_order,
        "source_order": None,
        "record_order": None,
        "source_id": None,
        "source_role": None,
        "raw_header": entry.raw_header,
        "id": entry.identifier,
        "description": entry.description,
        "sequence": entry.sequence,
        "kind": entry.kind,
        "contaminant_group": entry.contaminant_group,
        "sequence_hash": entry.sequence_hash,
        "entrapment_strategy": entry.entrapment_strategy,
        "decoy_strategy": entry.decoy_strategy,
        "decoy_source_order": entry.decoy_source_order,
        "decoy_source_id": entry.decoy_source_id,
        "decoy_source_kind": entry.decoy_source_kind,
    }


def _validate_frame(frame: pl.DataFrame, schema: pl.Schema, label: str) -> None:
    if frame.schema != schema:
        raise ValueError(f"invalid {label} schema: {frame.schema!r}")


def _optional_string(value: object) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: object) -> int | None:
    return None if value is None else int(str(value))
