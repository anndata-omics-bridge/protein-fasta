"""Independent biological-inventory to search-database workflow tests."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from protein_fasta import cli
from protein_fasta.analytics.hashing import file_checksum, sequence_hash
from protein_fasta.decoy_database import resolve_decoy_request, run_decoy_generation
from protein_fasta.inventory import (
    protein_inventory_frame,
    read_protein_inventory,
    read_search_inventory,
)
from protein_fasta.reading.parser import read_records
from protein_fasta.schema.decoy import DecoyRequestDocument


def _row(
    order: int,
    raw_header: str,
    sequence: str,
    kind: str,
    *,
    contaminant_group: str | None = None,
) -> dict[str, object]:
    identifier, separator, description = raw_header.partition(" ")
    has_source = kind in {"target", "contaminant", "entrapment"}
    return {
        "final_order": order,
        "source_order": order if has_source else None,
        "record_order": 0 if has_source else None,
        "source_id": identifier if has_source else None,
        "source_role": kind if has_source else None,
        "raw_header": raw_header,
        "id": identifier,
        "description": description if separator else None,
        "sequence": sequence,
        "kind": kind,
        "contaminant_group": contaminant_group,
        "sequence_hash": sequence_hash(sequence).hex(),
        "entrapment_strategy": None,
    }


def _inventory(path: Path) -> Path:
    protein_inventory_frame(
        [
            _row(0, "aa|demo|2026-08-28 biological", "CRAPCRAPCRAP", "sentinel"),
            _row(1, "sp|P1|ONE target", "MPEPTIDEK", "target"),
            _row(2, "aa|Cont_routine|routine contaminants", "MRECRAPCRAPCRAP", "sentinel"),
            _row(
                3,
                "sp|Cont_C1|C1 contaminant",
                "SAMPLEK",
                "contaminant",
                contaminant_group="routine",
            ),
        ]
    ).write_parquet(path)
    return path


def _request(output: Path) -> DecoyRequestDocument:
    return DecoyRequestDocument.model_validate(
        {
            "output_fasta": output,
            "strategy": {"type": "reverse"},
        }
    )


def test_reverse_strategy_rejects_irrelevant_seed() -> None:
    with pytest.raises(ValidationError, match="seed"):
        DecoyRequestDocument.model_validate(
            {
                "output_fasta": "search.fasta",
                "strategy": {"type": "reverse", "seed": 7},
            }
        )


def test_decoy_preserves_biological_order_and_links_every_generated_row(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "biological.parquet")
    checksum_before = file_checksum(inventory)
    effective = resolve_decoy_request(
        _request(Path("search.fasta")),
        request_base=tmp_path,
    )

    execution = run_decoy_generation(inventory, effective)

    assert file_checksum(inventory) == checksum_before
    records = list(read_records(execution.search_fasta_path))
    assert [record.raw_header for record in records[1:]] == [
        "sp|P1|ONE target",
        "aa|Cont_routine|routine contaminants",
        "sp|Cont_C1|C1 contaminant",
        "REV_sp|P1|ONE target",
        "REV_sp|Cont_C1|C1 contaminant",
    ]
    assert "; decoys reverse" in records[0].raw_header
    search = read_search_inventory(execution.search_inventory_path)
    assert search["kind"].to_list() == [
        "sentinel",
        "target",
        "sentinel",
        "contaminant",
        "decoy",
        "decoy",
    ]
    assert search["decoy_source_order"].drop_nulls().to_list() == [1, 3]
    assert search["decoy_source_id"].drop_nulls().to_list() == [
        "sp|P1|ONE",
        "sp|Cont_C1|C1",
    ]
    assert execution.document.counts.biological == 4
    assert execution.document.counts.decoy == 2
    assert execution.document.summary.n_sequences == 4
    assert execution.document.biological_inventory.checksum == checksum_before


def test_multiple_decoy_requests_reuse_one_unchanged_inventory(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "biological.parquet")
    first = run_decoy_generation(
        inventory,
        resolve_decoy_request(_request(Path("first.fasta")), request_base=tmp_path),
    )
    second = run_decoy_generation(
        inventory,
        resolve_decoy_request(_request(Path("second.fasta")), request_base=tmp_path),
    )

    assert (
        first.document.biological_inventory.checksum
        == second.document.biological_inventory.checksum
    )
    assert read_protein_inventory(inventory).height == 4


def test_decoy_rejects_an_inventory_that_already_contains_decoys(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "biological.parquet")
    frame = read_protein_inventory(inventory)
    frame = frame.vstack(
        protein_inventory_frame([_row(4, "REV_sp|P1|ONE target", "KEDITPEPM", "decoy")])
    )
    frame.write_parquet(inventory)

    with pytest.raises(ValueError, match="already contains decoy rows"):
        run_decoy_generation(
            inventory,
            resolve_decoy_request(_request(Path("search.fasta")), request_base=tmp_path),
        )


def test_cli_and_api_decoy_produce_identical_search_rows(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "biological.parquet")
    api = run_decoy_generation(
        inventory,
        resolve_decoy_request(_request(Path("api.fasta")), request_base=tmp_path),
    )
    request_path = tmp_path / "decoy.json"
    request_path.write_text(
        json.dumps(_request(Path("cli.fasta")).model_dump(mode="json")),
        encoding="utf-8",
    )

    cli.decoy(inventory, request=request_path)

    api_frame = read_search_inventory(api.search_inventory_path)
    cli_frame = read_search_inventory(tmp_path / "cli.fasta.search-inventory.parquet")
    assert cli_frame.schema == api_frame.schema
    assert cli_frame.rows() == api_frame.rows()


def test_decoy_direct_authors_request_and_replays(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "biological.parquet")
    output = tmp_path / "direct.fasta"

    cli.decoy(inventory, output=output, method="reverse")

    request_path = output.with_suffix(".fasta.request.json")
    assert json.loads(request_path.read_text(encoding="utf-8")) == {
        "decoy_prefix": "REV_",
        "output_fasta": "direct.fasta",
        "schema_version": "0.1",
        "strategy": {"type": "reverse"},
    }
    replay = tmp_path / "replay.fasta"
    cli.decoy(inventory, request=request_path, output=replay)
    assert read_search_inventory(output.with_suffix(".fasta.search-inventory.parquet")).rows() == (
        read_search_inventory(replay.with_suffix(".fasta.search-inventory.parquet")).rows()
    )
