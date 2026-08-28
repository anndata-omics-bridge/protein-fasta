"""Direct canonical-inventory indexing through the public API and CLI."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import pytest

from protein_fasta import cli
from protein_fasta.analytics.hashing import sequence_hash
from protein_fasta.database.models import BiologicalDatabase, ProteinInventoryEntry
from protein_fasta.inventory import biological_database_frame
from protein_fasta.registry.indexing import connect_registry, list_databases
from protein_fasta.registry_workflow import index_database_inventory
from protein_fasta.schema.registry import RegistryBackendDocument, RegistryDocument


def _write_inventory(tmp_path: Path) -> Path:
    path = tmp_path / "p1_db2_inventory_20260828.fasta.protein-inventory.parquet"
    entries = (
        ProteinInventoryEntry(
            final_order=0,
            source_order=None,
            record_order=None,
            source_id=None,
            source_role=None,
            raw_header="aa|p1_db2_inventory generated database",
            identifier="aa|p1_db2_inventory",
            description="generated database",
            sequence="CRAP",
            kind="sentinel",
            contaminant_group=None,
            sequence_hash=sequence_hash("CRAP").hex(),
            entrapment_strategy=None,
        ),
        ProteinInventoryEntry(
            final_order=1,
            source_order=0,
            record_order=0,
            source_id="human",
            source_role="target",
            raw_header="sp|P1|ONE first target",
            identifier="sp|P1|ONE",
            description="first target",
            sequence="PEPTIDE",
            kind="target",
            contaminant_group=None,
            sequence_hash=sequence_hash("PEPTIDE").hex(),
            entrapment_strategy=None,
        ),
        ProteinInventoryEntry(
            final_order=2,
            source_order=1,
            record_order=0,
            source_id="routine",
            source_role="contaminant",
            raw_header="sp|Cont_C1|C1 contaminant",
            identifier="sp|Cont_C1|C1",
            description="contaminant",
            sequence="SAMPLEK",
            kind="contaminant",
            contaminant_group="routine",
            sequence_hash=sequence_hash("SAMPLEK").hex(),
            entrapment_strategy=None,
        ),
    )
    biological_database_frame(BiologicalDatabase(entries)).write_parquet(path)
    return path


@pytest.mark.parametrize("backend", ["sqlite", "duckdb"])
def test_index_database_inventory_supports_both_backends(
    tmp_path: Path,
    backend: Literal["sqlite", "duckdb"],
) -> None:
    inventory = _write_inventory(tmp_path)
    suffix = ".sqlite3" if backend == "sqlite" else ".duckdb"
    registry_path = tmp_path / f"registry{suffix}"
    document = RegistryDocument(
        registry=RegistryBackendDocument(backend=backend),
    )

    record = index_database_inventory(inventory, registry_path, document)

    assert record.filename == "p1_db2_inventory_20260828.fasta"
    assert record.entry_count == 3
    assert record.target_count == 1
    assert record.contaminant_count == 1
    with connect_registry(registry_path, read_only=True) as connection:
        persisted = list_databases(connection)
    assert len(persisted) == 1
    assert persisted[0].target_content_fingerprint == record.target_content_fingerprint


def test_cli_and_api_inventory_index_have_identical_scientific_records(tmp_path: Path) -> None:
    inventory = _write_inventory(tmp_path)
    api_path = tmp_path / "api.sqlite3"
    cli_path = tmp_path / "cli.sqlite3"

    expected = index_database_inventory(inventory, api_path, RegistryDocument())
    cli.index_inventory(inventory, cli_path)
    with connect_registry(cli_path, read_only=True) as connection:
        actual = list_databases(connection)[0]

    assert actual.entry_count == expected.entry_count
    assert actual.kind_stats == expected.kind_stats
    assert actual.target_content_fingerprint == expected.target_content_fingerprint
