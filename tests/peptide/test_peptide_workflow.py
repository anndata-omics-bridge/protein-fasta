"""Protein-inventory peptide products are exact across execution backends."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import pytest

from protein_fasta import cli
from protein_fasta.analytics.hashing import sequence_hash
from protein_fasta.analytics_compile import make_digestion
from protein_fasta.inventory import search_inventory_frame
from protein_fasta.peptide.computation import (
    digest_partitions,
    partition_digest_inputs,
    peptide_database_from_partitions,
)
from protein_fasta.peptide.models import PeptideProtein
from protein_fasta.peptide_workflow import (
    read_peptide_comparisons,
    read_peptides,
    read_protein_peptide_map,
    resolve_peptide_build_request,
    resolve_peptide_comparison_request,
    run_peptide_build,
    run_peptide_comparison,
)
from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.peptide import (
    PeptideBuildRequestDocument,
    PeptideComparisonRequestDocument,
)


def _protein_row(
    order: int,
    identifier: str,
    sequence: str,
    kind: Literal["sentinel", "target", "contaminant", "entrapment"],
) -> dict[str, object]:
    return {
        "final_order": order,
        "source_order": order if kind != "sentinel" else None,
        "record_order": 0 if kind != "sentinel" else None,
        "source_id": "source" if kind != "sentinel" else None,
        "source_role": kind if kind != "sentinel" else None,
        "raw_header": identifier,
        "id": identifier,
        "description": None,
        "sequence": sequence,
        "kind": kind,
        "contaminant_group": "routine" if kind == "contaminant" else None,
        "sequence_hash": sequence_hash(sequence).hex(),
        "entrapment_strategy": "shuffled" if kind == "entrapment" else None,
        "decoy_strategy": None,
        "decoy_source_order": None,
        "decoy_source_id": None,
        "decoy_source_kind": None,
    }


def _inventory(path: Path) -> Path:
    rows = [
        _protein_row(0, "aa|demo|metadata", "CRAP", "sentinel"),
        _protein_row(1, "sp|P1|ONE", "MPEPTIDEKAAAK", "target"),
        _protein_row(2, "sp|C1|CONT", "MPEPTIDEK", "contaminant"),
        _protein_row(3, "sp|E1|ENTRAP", "MENTRAPK", "entrapment"),
    ]
    rows.append(
        {
            "final_order": 4,
            "source_order": None,
            "record_order": None,
            "source_id": None,
            "source_role": None,
            "raw_header": "REV_sp|P1|ONE",
            "id": "REV_sp|P1|ONE",
            "description": None,
            "sequence": "KAAAEDITPEPM",
            "kind": "decoy",
            "contaminant_group": None,
            "sequence_hash": sequence_hash("KAAAEDITPEPM").hex(),
            "entrapment_strategy": None,
            "decoy_strategy": "reverse",
            "decoy_source_order": 1,
            "decoy_source_id": "sp|P1|ONE",
            "decoy_source_kind": "target",
        }
    )
    search_inventory_frame(rows).write_parquet(path)
    return path


def _request(
    prefix: str,
    backend: Literal["memory", "sqlite", "duckdb"],
) -> PeptideBuildRequestDocument:
    return PeptideBuildRequestDocument.model_validate(
        {
            "peptides_parquet": f"{prefix}-peptides.parquet",
            "mapping_parquet": f"{prefix}-mapping.parquet",
            "peptide_fasta": f"{prefix}-peptides.fasta",
            "digestion": {
                "enzyme": "trypsin",
                "min_length": 3,
                "max_length": 30,
                "missed_cleavages": 1,
            },
            "execution": {"type": backend, "workers": 1, "partition_size": 2},
        }
    )


def test_memory_sqlite_and_duckdb_produce_identical_peptide_artifacts(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "search.parquet")
    executions = [
        run_peptide_build(
            inventory,
            resolve_peptide_build_request(_request(backend, backend), request_base=tmp_path),
        )
        for backend in ("memory", "sqlite", "duckdb")
    ]

    peptide_rows = [read_peptides(execution.peptides_path).rows() for execution in executions]
    mapping_rows = [
        read_protein_peptide_map(execution.mapping_path).rows() for execution in executions
    ]
    assert peptide_rows[0] == peptide_rows[1] == peptide_rows[2]
    assert mapping_rows[0] == mapping_rows[1] == mapping_rows[2]
    assert executions[0].document.counts.input_proteins == 4
    assert executions[0].document.counts.entrapment_peptides > 0
    assert executions[0].document.counts.decoy_peptides > 0


def test_peptides_cli_matches_api_and_comparison_is_exact(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "search.parquet")
    api = run_peptide_build(
        inventory,
        resolve_peptide_build_request(_request("api", "memory"), request_base=tmp_path),
    )
    request_path = tmp_path / "peptides.json"
    request_path.write_text(
        json.dumps(_request("cli", "memory").model_dump(mode="json")),
        encoding="utf-8",
    )

    cli.peptides(inventory, request=request_path)

    cli_peptides = tmp_path / "cli-peptides.parquet"
    assert read_peptides(api.peptides_path).rows() == read_peptides(cli_peptides).rows()
    comparison = run_peptide_comparison(
        api.peptides_path,
        cli_peptides,
        resolve_peptide_comparison_request(
            PeptideComparisonRequestDocument(output_parquet=Path("comparison.parquet")),
            request_base=tmp_path,
        ),
    )
    rows = read_peptide_comparisons(comparison.comparison_path)
    assert rows["jaccard"].to_list() == [1.0] * 5


def test_peptide_commands_author_requests_and_replay_with_new_outputs(tmp_path: Path) -> None:
    inventory = _inventory(tmp_path / "search.parquet")
    direct_dir = tmp_path / "direct"

    cli.peptides(inventory, output=direct_dir, minimum=3, maximum=30, missed=1)

    request_path = direct_dir / "peptides.request.json"
    payload = json.loads(request_path.read_text(encoding="utf-8"))
    assert payload == {
        "digestion": {
            "enzyme": "trypsin",
            "max_length": 30,
            "min_length": 3,
            "missed_cleavages": 1,
        },
        "execution": {"partition_size": 500, "type": "memory", "workers": 1},
        "mapping_parquet": "protein-peptide-map.parquet",
        "peptide_fasta": "peptides.fasta",
        "peptides_parquet": "peptides.parquet",
        "schema_version": "0.1",
    }
    replay_dir = tmp_path / "replay"
    cli.peptides(inventory, request=request_path, output=replay_dir)
    assert (
        read_peptides(direct_dir / "peptides.parquet").rows()
        == read_peptides(replay_dir / "peptides.parquet").rows()
    )

    saved_dir = tmp_path / "saved"
    saved_request = tmp_path / "requests" / "peptides.json"
    cli.peptides(inventory, output=saved_dir, save=saved_request)
    saved_payload = json.loads(saved_request.read_text(encoding="utf-8"))
    assert saved_payload["peptides_parquet"] == "../saved/peptides.parquet"
    assert (saved_dir / "peptides.parquet").is_file()

    comparison = tmp_path / "direct-comparison.parquet"
    cli.pepcompare(
        direct_dir / "peptides.parquet",
        replay_dir / "peptides.parquet",
        output=comparison,
    )
    comparison_request = comparison.with_suffix(".parquet.request.json")
    assert json.loads(comparison_request.read_text(encoding="utf-8")) == {
        "output_parquet": "direct-comparison.parquet",
        "schema_version": "0.1",
    }
    replay_comparison = tmp_path / "replay-comparison.parquet"
    cli.pepcompare(
        direct_dir / "peptides.parquet",
        replay_dir / "peptides.parquet",
        request=comparison_request,
        output=replay_comparison,
    )
    assert (
        read_peptide_comparisons(comparison).rows()
        == read_peptide_comparisons(replay_comparison).rows()
    )


def test_repeated_protein_identifier_produces_one_canonical_mapping() -> None:
    proteins = (
        PeptideProtein(0, "P1", "target", "AK"),
        PeptideProtein(1, "P1", "target", "AK"),
    )
    tasks = partition_digest_inputs(
        proteins,
        make_digestion(DigestionDocument(min_length=1, max_length=20)),
        partition_size=1,
    )

    database = peptide_database_from_partitions(digest_partitions(tasks, workers=1))

    assert len(database.mappings) == 1
    assert database.peptides[0].mapping_count == 1
    assert database.peptides[0].protein_count == 1


def test_conflicting_repeated_protein_mapping_is_rejected() -> None:
    proteins = (
        PeptideProtein(0, "P1", "target", "AK"),
        PeptideProtein(1, "P1", "contaminant", "AK"),
    )
    tasks = partition_digest_inputs(
        proteins,
        make_digestion(DigestionDocument(min_length=1, max_length=20)),
        partition_size=1,
    )

    with pytest.raises(ValueError, match="conflicting protein kind"):
        peptide_database_from_partitions(digest_partitions(tasks, workers=1))
