"""Decoy-method diagnostics reuse one unchanged biological inventory."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from protein_fasta import cli
from protein_fasta.analytics.hashing import file_checksum, sequence_hash
from protein_fasta.decoy_report import (
    read_decoy_report,
    resolve_decoy_report_request,
    run_decoy_report,
)
from protein_fasta.inventory import protein_inventory_frame
from protein_fasta.schema.decoy_report import DecoyReportRequestDocument

pytestmark = pytest.mark.generation


def _inventory(path: Path) -> Path:
    rows: list[dict[str, object]] = []
    for order, (identifier, sequence, kind) in enumerate(
        (
            ("aa|demo|metadata", "CRAPCRAP", "sentinel"),
            ("sp|P1|ONE", "MPEPTIDEKAAALEVSSK", "target"),
            ("sp|P2|TWO", "MSAMPLERGGGYTFDK", "target"),
        )
    ):
        rows.append(
            {
                "final_order": order,
                "source_order": None if kind == "sentinel" else 0,
                "record_order": None if kind == "sentinel" else order - 1,
                "source_id": None if kind == "sentinel" else "targets",
                "source_role": None if kind == "sentinel" else "target",
                "raw_header": identifier,
                "id": identifier,
                "description": None,
                "sequence": sequence,
                "kind": kind,
                "contaminant_group": None,
                "sequence_hash": sequence_hash(sequence).hex(),
                "entrapment_strategy": None,
            }
        )
    protein_inventory_frame(rows).write_parquet(path)
    return path


def _request(output: str) -> DecoyReportRequestDocument:
    return DecoyReportRequestDocument.model_validate(
        {
            "output_parquet": output,
            "digestion": {
                "enzyme": "trypsin",
                "min_length": 3,
                "max_length": 30,
                "missed_cleavages": 1,
            },
            "strategies": [
                {"type": "reverse"},
                {"type": "shuffle", "seed": 7},
                {
                    "type": "decoypyrat",
                    "seed": 7,
                    "digestion": {
                        "enzyme": "trypsin",
                        "min_length": 3,
                        "max_length": 30,
                        "missed_cleavages": 1,
                    },
                },
            ],
        }
    )


def test_all_methods_compare_against_one_unchanged_inventory(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "biological.parquet")
    checksum = file_checksum(inventory)

    execution = run_decoy_report(
        inventory,
        resolve_decoy_report_request(_request("report.parquet"), request_base=tmp_path),
    )

    assert file_checksum(inventory) == checksum
    frame = read_decoy_report(execution.comparison_path)
    assert frame["method"].to_list() == [
        "biological",
        "reverse",
        "shuffle",
        "decoypyrat",
    ]
    assert execution.document.biological_inventory.checksum == checksum
    assert execution.document.target_proteins == 2


def test_decoy_report_cli_matches_api(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "biological.parquet")
    api = run_decoy_report(
        inventory,
        resolve_decoy_report_request(_request("api.parquet"), request_base=tmp_path),
    )
    request_path = tmp_path / "report.json"
    request_path.write_text(
        json.dumps(_request("cli.parquet").model_dump(mode="json")),
        encoding="utf-8",
    )

    cli.decoy_report(inventory, request_path)

    assert (
        read_decoy_report(api.comparison_path).rows()
        == read_decoy_report(tmp_path / "cli.parquet").rows()
    )
