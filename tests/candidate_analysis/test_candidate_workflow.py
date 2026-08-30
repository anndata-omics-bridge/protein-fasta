"""Canonical-inventory candidate review is read-only and backend-equivalent."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import pytest

from protein_fasta import cli
from protein_fasta.analytics.hashing import file_checksum, sequence_hash
from protein_fasta.candidate_analysis import (
    read_candidate_comparisons,
    resolve_candidate_request,
    run_candidate_analysis,
)
from protein_fasta.inventory import protein_inventory_frame
from protein_fasta.registry.indexing import rebuild_registry
from protein_fasta.registry_workflow import make_registry_settings
from protein_fasta.schema.candidate import CandidateRequestDocument
from protein_fasta.schema.registry import RegistryBackendDocument, RegistryDocument


def _inventory(path: Path) -> Path:
    protein_inventory_frame(
        [
            {
                "final_order": 0,
                "source_order": 0,
                "record_order": 0,
                "source_id": "selected",
                "source_role": "target",
                "raw_header": "sp|P1|ONE candidate",
                "id": "sp|P1|ONE",
                "description": "candidate",
                "sequence": "MPEPTIDEK",
                "kind": "target",
                "contaminant_group": None,
                "sequence_hash": sequence_hash("MPEPTIDEK").hex(),
                "entrapment_strategy": None,
            }
        ]
    ).write_parquet(path)
    return path


def _registry(
    tmp_path: Path,
    backend: Literal["sqlite", "duckdb"],
) -> tuple[Path, RegistryDocument]:
    fasta_root = tmp_path / "registered"
    fasta_root.mkdir()
    (fasta_root / "registered_20260828.fasta").write_text(
        ">sp|P1|ONE registered\nMPEPTIDEK\n",
        encoding="utf-8",
    )
    suffix = ".sqlite3" if backend == "sqlite" else ".duckdb"
    path = tmp_path / f"registry{suffix}"
    document = RegistryDocument(registry=RegistryBackendDocument(backend=backend))
    rebuild_registry(
        fasta_root,
        path,
        make_registry_settings(document, fasta_root=fasta_root, registry_path=path),
    )
    return path, document


@pytest.mark.parametrize("backend", ["sqlite", "duckdb"])
def test_candidate_review_reads_inventory_without_mutating_registry(
    tmp_path: Path,
    backend: Literal["sqlite", "duckdb"],
) -> None:
    inventory = _inventory(tmp_path / "candidate.parquet")
    registry, document = _registry(tmp_path, backend)
    checksum_before = file_checksum(registry)
    effective = resolve_candidate_request(
        CandidateRequestDocument(output_parquet=Path("comparison.parquet")),
        request_base=tmp_path,
    )

    execution = run_candidate_analysis(inventory, registry, effective, document)

    assert file_checksum(registry) == checksum_before
    comparison = read_candidate_comparisons(execution.comparison_path)
    assert comparison["kind"].to_list() == ["target", "contaminant"]
    assert comparison["relationship"].to_list()[0] == "exact_content"
    assert execution.document.counts.checked_databases == 1
    assert execution.analysis.neighbourhood.relative_paths == (
        "candidate.parquet [candidate]",
        "registered_20260828.fasta",
    )


def test_candidate_cli_and_api_write_identical_comparison_rows(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "candidate.parquet")
    registry, document = _registry(tmp_path, "sqlite")
    api = run_candidate_analysis(
        inventory,
        registry,
        resolve_candidate_request(
            CandidateRequestDocument(output_parquet=Path("api.parquet")),
            request_base=tmp_path,
        ),
        document,
    )
    request_path = tmp_path / "candidate.json"
    request_path.write_text(
        json.dumps(
            CandidateRequestDocument(output_parquet=Path("cli.parquet")).model_dump(mode="json")
        ),
        encoding="utf-8",
    )
    config_path = tmp_path / "registry.json"
    config_path.write_text(json.dumps(document.model_dump(mode="json")), encoding="utf-8")

    cli.candidate(inventory, registry, request=request_path, config=config_path)

    assert (
        read_candidate_comparisons(api.comparison_path).rows()
        == read_candidate_comparisons(tmp_path / "cli.parquet").rows()
    )


def test_candidate_direct_authors_request_and_replays(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "candidate.parquet")
    registry, document = _registry(tmp_path, "sqlite")
    config_path = tmp_path / "registry.json"
    config_path.write_text(json.dumps(document.model_dump(mode="json")), encoding="utf-8")
    output = tmp_path / "direct.parquet"

    cli.candidate(inventory, registry, output=output, config=config_path, limit=10)

    request_path = output.with_suffix(".parquet.request.json")
    assert json.loads(request_path.read_text(encoding="utf-8")) == {
        "clustering_metric": "target_ids",
        "neighbour_limit": 10,
        "output_parquet": "direct.parquet",
        "overlap_threshold": 0.99,
        "schema_version": "0.1",
    }
    replay = tmp_path / "replay.parquet"
    cli.candidate(
        inventory,
        registry,
        request=request_path,
        output=replay,
        config=config_path,
    )
    assert read_candidate_comparisons(output).rows() == read_candidate_comparisons(replay).rows()
