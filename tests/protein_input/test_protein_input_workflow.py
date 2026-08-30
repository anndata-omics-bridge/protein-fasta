"""Canonical protein-input preparation through API and CLI."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from protein_fasta import cli
from protein_fasta.analytics.hashing import file_checksum
from protein_fasta.database.models import (
    BiologicalDatabase,
    BiologicalKind,
    DecoyInventoryEntry,
    ProteinInventoryEntry,
    SearchDatabase,
)
from protein_fasta.inventory import (
    biological_database_frame,
    read_protein_input,
    search_database_frame,
)
from protein_fasta.protein_input import (
    ContaminantProteinBlock,
    ProteinSourceEntry,
    TargetProteinSource,
    prepare_protein_input_frame,
    resolve_derived_protein_input_request,
    resolve_protein_input_request,
    run_derived_protein_input_preparation,
    run_protein_input_preparation,
)
from protein_fasta.schema.protein_input import (
    DerivedProteinInputRequestDocument,
    ProteinInputRequestDocument,
)


def _request(output: Path) -> ProteinInputRequestDocument:
    return ProteinInputRequestDocument.model_validate(
        {
            "sources": [
                {"type": "target", "source_id": "human", "path": "targets.fasta"},
                {
                    "type": "contaminant",
                    "source_id": "routine",
                    "path": "contaminants.fasta",
                    "block_name": "routine",
                    "block_description": "routine contaminants",
                },
            ],
            "output_parquet": output,
        }
    )


def _write_sources(directory: Path) -> None:
    directory.joinpath("targets.fasta").write_text(
        ">sp|P1|ONE target one\nak*\n>sp|P2|TWO target two\nPEPTIDE\n",
        encoding="utf-8",
    )
    directory.joinpath("contaminants.fasta").write_text(
        ">sp|Cont_C1|C1 contaminant\nSAMPLEK\n",
        encoding="utf-8",
    )


def test_source_variants_reject_contaminant_without_a_block_name() -> None:
    with pytest.raises(ValidationError, match="block_name"):
        ProteinInputRequestDocument.model_validate(
            {
                "sources": [
                    {"type": "contaminant", "source_id": "routine", "path": "source.fasta"}
                ],
                "output_parquet": "input.parquet",
            }
        )


def test_prepare_preserves_order_roles_normalization_and_exact_evidence(tmp_path: Path) -> None:
    _write_sources(tmp_path)
    effective = resolve_protein_input_request(
        _request(Path("protein-input.parquet")),
        request_base=tmp_path,
    )

    execution = run_protein_input_preparation(effective)

    frame = read_protein_input(execution.protein_input_path)
    assert list(
        zip(frame["source_id"].to_list(), frame["record_order"].to_list(), strict=True)
    ) == [
        ("human", 0),
        ("human", 1),
        ("routine", 0),
    ]
    assert frame["role"].to_list() == ["target", "target", "contaminant"]
    assert frame["block_name"].to_list() == [None, None, "routine"]
    assert frame["sequence"].to_list() == ["AK", "PEPTIDE", "SAMPLEK"]
    assert execution.document.normalization.upper_cased == 1
    assert execution.document.normalization.terminal_stops_stripped == 1
    assert execution.document.protein_input.checksum == file_checksum(execution.protein_input_path)
    assert [source.artifact.row_count for source in execution.document.sources] == [2, 1]
    assert execution.effective_request_path.is_file()
    assert execution.result_path.is_file()


def test_in_memory_preparation_matches_the_file_adapter(tmp_path: Path) -> None:
    _write_sources(tmp_path)
    file_execution = run_protein_input_preparation(
        resolve_protein_input_request(
            _request(Path("protein-input.parquet")),
            request_base=tmp_path,
        )
    )

    prepared = prepare_protein_input_frame(
        (
            TargetProteinSource(
                "human",
                (
                    ProteinSourceEntry("sp|P1|ONE target one", "ak*"),
                    ProteinSourceEntry("sp|P2|TWO target two", "PEPTIDE"),
                ),
            ),
        ),
        (
            ContaminantProteinBlock(
                "routine",
                "routine",
                "routine contaminants",
                (ProteinSourceEntry("sp|Cont_C1|C1 contaminant", "SAMPLEK"),),
            ),
        ),
    )

    assert prepared.frame.equals(file_execution.frame)
    assert prepared.source_row_counts == (2, 1)
    assert prepared.upper_cased == 1
    assert prepared.terminal_stops_stripped == 1


def test_failed_prepare_keeps_effective_request_without_success_artifacts(tmp_path: Path) -> None:
    effective = resolve_protein_input_request(
        ProteinInputRequestDocument.model_validate(
            {
                "sources": [{"type": "target", "source_id": "missing", "path": "missing.fasta"}],
                "output_parquet": "protein-input.parquet",
            }
        ),
        request_base=tmp_path,
    )

    with pytest.raises(ValueError, match=r"missing\.fasta"):
        run_protein_input_preparation(effective)

    assert effective.output_parquet.with_suffix(".parquet.effective.json").is_file()
    assert not effective.output_parquet.exists()
    assert not effective.output_parquet.with_suffix(".parquet.result.json").exists()


def test_cli_and_api_prepare_equivalent_frames(tmp_path: Path) -> None:
    _write_sources(tmp_path)
    api_request = _request(Path("api.parquet"))
    api_execution = run_protein_input_preparation(
        resolve_protein_input_request(api_request, request_base=tmp_path)
    )
    cli_request = _request(Path("cli.parquet"))
    request_path = tmp_path / "prepare.json"
    request_path.write_text(
        json.dumps(cli_request.model_dump(mode="json")),
        encoding="utf-8",
    )

    cli.prepare(request=request_path)

    api_frame = read_protein_input(api_execution.protein_input_path)
    cli_frame = read_protein_input(tmp_path / "cli.parquet")
    assert cli_frame.schema == api_frame.schema
    assert cli_frame.rows() == api_frame.rows()


def test_prepare_direct_authors_request_and_replays_with_output_override(tmp_path: Path) -> None:
    _write_sources(tmp_path)
    output = tmp_path / "direct.parquet"

    cli.prepare(tmp_path / "targets.fasta", output, id="human")

    request_path = output.with_suffix(".parquet.request.json")
    assert json.loads(request_path.read_text(encoding="utf-8")) == {
        "output_parquet": "direct.parquet",
        "schema_version": "0.1",
        "sources": [{"path": "targets.fasta", "source_id": "human", "type": "target"}],
    }
    replay = tmp_path / "replay.parquet"
    cli.prepare(request=request_path, output=replay)
    assert read_protein_input(output).rows() == read_protein_input(replay).rows()

    with pytest.raises(FileExistsError):
        cli.prepare(tmp_path / "targets.fasta", output, id="human")
    with pytest.raises(ValueError, match="cannot be combined with --request"):
        cli.prepare(tmp_path / "targets.fasta", request=request_path, id="human")


def _protein(
    order: int,
    identifier: str,
    kind: BiologicalKind,
    *,
    group: str | None = None,
) -> ProteinInventoryEntry:
    return ProteinInventoryEntry(
        final_order=order,
        source_order=None,
        record_order=None,
        source_id=None,
        source_role=None,
        raw_header=identifier,
        identifier=identifier,
        description=None,
        sequence="PEPTIDE",
        kind=kind,
        contaminant_group=group,
        sequence_hash="hash",
        entrapment_strategy="shuffled" if kind == "entrapment" else None,
    )


def _write_derived_inventories(tmp_path: Path) -> tuple[Path, Path]:
    source_path = tmp_path / "source.search-inventory.parquet"
    source_entries = (
        _protein(0, "aa|source", "sentinel"),
        _protein(1, "P1", "target"),
        _protein(2, "Cont_C1", "contaminant", group="routine"),
        _protein(3, "ENTRAP_P1", "entrapment"),
        DecoyInventoryEntry(
            final_order=4,
            raw_header="rev_P1",
            identifier="rev_P1",
            description=None,
            sequence="EDITPEP",
            contaminant_group=None,
            sequence_hash="decoy-hash",
            entrapment_strategy=None,
            decoy_strategy="reverse",
            decoy_source_order=1,
            decoy_source_id="P1",
            decoy_source_kind="target",
        ),
    )
    search_database_frame(SearchDatabase(source_entries)).write_parquet(source_path)

    foreign_path = tmp_path / "foreign.protein-inventory.parquet"
    biological_database_frame(
        BiologicalDatabase(
            (
                _protein(0, "aa|foreign", "sentinel"),
                _protein(1, "F1", "target"),
                _protein(2, "Cont_F1", "contaminant", group="foreign-routine"),
                _protein(3, "ENTRAP_F1", "entrapment"),
            )
        )
    ).write_parquet(foreign_path)
    return source_path, foreign_path


def _derived_request(output: str) -> DerivedProteinInputRequestDocument:
    return DerivedProteinInputRequestDocument(
        source_inventory=Path("source.search-inventory.parquet"),
        source_id="source-db",
        foreign_inventory=Path("foreign.protein-inventory.parquet"),
        foreign_source_id="foreign-db",
        output_parquet=Path(output),
    )


def test_derive_input_keeps_only_original_targets_and_contaminants(tmp_path: Path) -> None:
    _write_derived_inventories(tmp_path)
    execution = run_derived_protein_input_preparation(
        resolve_derived_protein_input_request(
            _derived_request("derived.parquet"),
            request_base=tmp_path,
        )
    )

    frame = read_protein_input(execution.protein_input_path)
    assert frame["id"].to_list() == ["P1", "Cont_C1", "F1", "Cont_F1"]
    assert frame["role"].to_list() == ["target", "contaminant", "foreign", "foreign"]
    assert frame["block_name"].to_list() == [None, "routine", None, None]
    assert execution.document.counts.model_dump() == {
        "target": 1,
        "contaminant": 1,
        "foreign": 2,
        "skipped_sentinel": 2,
        "skipped_entrapment": 2,
        "skipped_decoy": 1,
    }
    assert [source.artifact.schema_name for source in execution.document.sources] == [
        "search-inventory",
        "protein-inventory",
    ]


def test_cli_and_api_derived_inputs_are_equivalent(tmp_path: Path) -> None:
    _write_derived_inventories(tmp_path)
    api_execution = run_derived_protein_input_preparation(
        resolve_derived_protein_input_request(
            _derived_request("api-derived.parquet"),
            request_base=tmp_path,
        )
    )
    request_path = tmp_path / "derive-input.json"
    request_path.write_text(
        json.dumps(_derived_request("cli-derived.parquet").model_dump(mode="json")),
        encoding="utf-8",
    )

    cli.derive_input(request=request_path)

    assert (
        read_protein_input(tmp_path / "cli-derived.parquet").rows()
        == read_protein_input(api_execution.protein_input_path).rows()
    )


def test_derive_input_direct_authors_request_and_replays(tmp_path: Path) -> None:
    source, _ = _write_derived_inventories(tmp_path)
    output = tmp_path / "direct-derived.parquet"

    cli.derive_input(source, output, id="source-db")

    request_path = output.with_suffix(".parquet.request.json")
    assert json.loads(request_path.read_text(encoding="utf-8")) == {
        "output_parquet": "direct-derived.parquet",
        "schema_version": "0.1",
        "source_id": "source-db",
        "source_inventory": "source.search-inventory.parquet",
    }
    replay = tmp_path / "replay-derived.parquet"
    cli.derive_input(request=request_path, output=replay)
    assert read_protein_input(output).rows() == read_protein_input(replay).rows()
