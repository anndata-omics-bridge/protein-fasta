"""Shared database assembly keeps construction separate from registry analysis."""

from __future__ import annotations

import datetime
from pathlib import Path

import polars as pl
import pytest

import protein_fasta.database_build as database_build
from protein_fasta.database.entrapment import EntrapmentStrategy
from protein_fasta.database_build import (
    BiologicalBuildResult,
    ConflictingIdError,
    ContaminantBlock,
    Entry,
    build_database,
    run_database_build,
    run_database_build_from_frame,
)
from protein_fasta.database_compile import (
    make_database_metadata,
    make_database_naming,
    make_entrapment_generation,
)
from protein_fasta.inventory import (
    PROTEIN_INPUT_SCHEMA,
    read_entrapment_pairs,
    read_protein_input,
)
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import read_records
from protein_fasta.registry.rules import load_registry_diagnostics
from protein_fasta.schema.build import (
    EffectiveDatabaseBuildDocument,
    ForeignSpeciesEntrapmentDocument,
    MetadataDocument,
    NamingDocument,
    ShuffledEntrapmentDocument,
)


def _build(
    tmp_path: Path,
    *,
    targets: tuple[Entry, ...] = (("sp|P1|ONE Protein one", "ak*"),),
    contaminant_blocks: tuple[ContaminantBlock, ...] = (),
) -> BiologicalBuildResult:
    return build_database(
        targets=targets,
        name_fields={"description": "demo"},
        naming=make_database_naming(NamingDocument(default_dbname="derived")),
        metadata=make_database_metadata(MetadataDocument()),
        diagnostics=load_registry_diagnostics(),
        output_dir=tmp_path,
        date=datetime.date(2026, 8, 27),
        template="derived",
        contaminant_blocks=contaminant_blocks,
    )


def _protein_input(path: Path, *, contaminant: bool = False) -> Path:
    rows: list[dict[str, object]] = [
        {
            "source_order": 0,
            "record_order": 0,
            "source_id": "targets",
            "role": "target",
            "block_name": None,
            "block_description": None,
            "raw_header": "sp|P1|ONE target",
            "id": "sp|P1|ONE",
            "description": "target",
            "sequence": "PEPTIDEK",
            "upper_cased": False,
            "terminal_stop_stripped": False,
        }
    ]
    if contaminant:
        rows.append(
            {
                "source_order": 1,
                "record_order": 0,
                "source_id": "routine",
                "role": "contaminant",
                "block_name": "routine",
                "block_description": "routine contaminants",
                "raw_header": "sp|Cont_C1|C1 contaminant",
                "id": "sp|Cont_C1|C1",
                "description": "contaminant",
                "sequence": "SAMPLEK",
                "upper_cased": False,
                "terminal_stop_stripped": False,
            }
        )
    pl.DataFrame(rows, schema=PROTEIN_INPUT_SCHEMA).write_parquet(path)
    return path


def test_effective_request_is_written_before_source_reading(tmp_path: Path) -> None:
    effective = EffectiveDatabaseBuildDocument(
        output_dir=tmp_path / "out",
        date=datetime.date(2026, 8, 28),
        name_fields={"description": "failed"},
        template="derived",
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
    )

    with pytest.raises((FileNotFoundError, OSError), match=r"missing\.parquet"):
        run_database_build(tmp_path / "missing.parquet", effective)

    effective_path = tmp_path / "out" / "failed_20260828.fasta.effective.json"
    assert effective_path.is_file()


def test_workflow_inventory_preserves_final_kinds_and_generation_evidence(tmp_path: Path) -> None:
    protein_input = _protein_input(tmp_path / "protein-input.parquet", contaminant=True)
    effective = EffectiveDatabaseBuildDocument(
        output_dir=tmp_path / "out",
        date=datetime.date(2026, 8, 28),
        name_fields={"description": "inventory"},
        template="derived",
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
    )

    execution = run_database_build(protein_input, effective)

    inventory = pl.read_parquet(effective.output_dir / execution.document.protein_inventory.path)
    assert inventory["kind"].to_list() == [
        "sentinel",
        "target",
        "sentinel",
        "contaminant",
    ]
    assert inventory["contaminant_group"].to_list()[3] == "routine"
    assert inventory["entrapment_strategy"].null_count() == inventory.height
    assert execution.document.counts.total == 4
    assert not hasattr(execution.document.counts, "decoy")


def test_frame_build_persists_replay_input_and_matches_artifact_build(tmp_path: Path) -> None:
    input_path = _protein_input(tmp_path / "protein-input.parquet", contaminant=True)
    frame = read_protein_input(input_path)
    path_effective = EffectiveDatabaseBuildDocument(
        output_dir=tmp_path / "path",
        date=datetime.date(2026, 8, 28),
        name_fields={"description": "parity"},
        template="derived",
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
    )
    frame_effective = path_effective.model_copy(update={"output_dir": tmp_path / "frame"})

    artifact_execution = run_database_build(input_path, path_effective)
    frame_execution = run_database_build_from_frame(frame, frame_effective)

    assert read_protein_input(frame_execution.protein_input_path).equals(frame)
    assert frame_execution.result.path.read_bytes() == artifact_execution.result.path.read_bytes()
    assert frame_execution.database == artifact_execution.database


def test_frame_build_rejects_noncanonical_schema_before_writing(tmp_path: Path) -> None:
    effective = EffectiveDatabaseBuildDocument(
        output_dir=tmp_path / "out",
        date=datetime.date(2026, 8, 28),
        name_fields={"description": "invalid"},
        template="derived",
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
    )

    with pytest.raises(ValueError, match="invalid protein-input schema"):
        run_database_build_from_frame(pl.DataFrame({"id": ["P1"]}), effective)

    assert not effective.output_dir.exists()


def test_result_publication_failure_rolls_back_biological_artifacts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protein_input = _protein_input(tmp_path / "protein-input.parquet")
    effective = EffectiveDatabaseBuildDocument(
        output_dir=tmp_path / "out",
        date=datetime.date(2026, 8, 28),
        name_fields={"description": "atomic"},
        template="derived",
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
    )
    real_write_json = database_build.write_json_atomic

    def fail_result(
        path: Path,
        payload: object,
        /,
        *,
        replace_existing: bool,
    ) -> None:
        if path.name.endswith(".result.json"):
            raise OSError("simulated result publication failure")
        real_write_json(path, payload, replace_existing=replace_existing)

    monkeypatch.setattr(database_build, "write_json_atomic", fail_result)

    with pytest.raises(OSError, match="simulated result publication failure"):
        run_database_build(protein_input, effective)

    prefix = effective.output_dir / "atomic_20260828.fasta"
    assert prefix.with_suffix(".fasta.effective.json").is_file()
    assert not prefix.exists()
    assert not prefix.with_suffix(".fasta.protein-inventory.parquet").exists()
    assert not prefix.with_suffix(".fasta.result.json").exists()


def test_build_writes_metadata_targets_and_contaminant_markers(tmp_path: Path) -> None:
    result = _build(
        tmp_path,
        contaminant_blocks=(
            ContaminantBlock("crap", "routine contaminants", (("sp|Cont_C1|C1", "PEPTIDE"),)),
        ),
    )

    records = list(read_records(result.path))
    assert records[0].raw_header.startswith("aa|demo|")
    assert [record.raw_header for record in records[1:]] == [
        "sp|P1|ONE Protein one",
        "aa|Cont_crap|routine contaminants",
        "sp|Cont_C1|C1",
    ]
    assert records[1].sequence == "AK"
    assert result.n_target == 1
    assert result.n_contaminant == 1
    assert result.upper_cased_entries == 1
    assert result.stop_stripped_entries == 1


def test_build_rejects_conflicting_ids_after_normalization(tmp_path: Path) -> None:
    with pytest.raises(ConflictingIdError, match="different sequences"):
        _build(
            tmp_path,
            targets=(("P1 first", "AAAA"), ("P1 second", "BBBB")),
        )


@pytest.mark.generation
def test_entrapment_joins_the_biological_database_without_decoys(tmp_path: Path) -> None:
    diagnostics = load_registry_diagnostics()
    result = build_database(
        targets=(
            ("sp|P1|ONE first", "MPEPTIDEKAAALEVSSKMWNQTR"),
            ("sp|P2|TWO second", "MSAMPLERGGGYTFDKQVNTLPWR"),
        ),
        name_fields={"description": "entrapment"},
        naming=make_database_naming(NamingDocument(default_dbname="derived")),
        metadata=make_database_metadata(MetadataDocument()),
        diagnostics=diagnostics,
        output_dir=tmp_path,
        date=datetime.date(2026, 8, 27),
        template="derived",
        entrapment_generation=make_entrapment_generation(
            ShuffledEntrapmentDocument(fold=1, seed=7)
        ),
    )

    identifiers = [parse_header(record.raw_header).id for record in read_records(result.path)]
    entrapment_ids = {
        identifier
        for identifier in identifiers
        if identifier.endswith("_p_target") and not identifier.startswith(diagnostics.decoy_prefix)
    }

    assert result.n_entrapment == 2
    assert len(entrapment_ids) == result.n_entrapment
    assert not {f"REV_{identifier}" for identifier in entrapment_ids} & set(identifiers)
    assert result.entrapment_pairs_path is not None
    assert result.entrapment_pairs_path.is_file()
    assert read_entrapment_pairs(result.entrapment_pairs_path).height > 0
    assert result.entrapment is not None
    assert result.entrapment.strategy == EntrapmentStrategy.SHUFFLED
    assert result.entrapment.requested_fold == result.entrapment.achieved_fold == 1


@pytest.mark.generation
def test_foreign_species_entrapment_does_not_claim_peptide_pairs(tmp_path: Path) -> None:
    result = build_database(
        targets=(("sp|P1|ONE first", "MPEPTIDEKTESTSEQUENCEK"),),
        name_fields={"description": "foreign"},
        naming=make_database_naming(NamingDocument(default_dbname="derived")),
        metadata=make_database_metadata(MetadataDocument()),
        diagnostics=load_registry_diagnostics(),
        output_dir=tmp_path,
        date=datetime.date(2026, 8, 27),
        template="derived",
        entrapment_generation=make_entrapment_generation(ForeignSpeciesEntrapmentDocument(fold=1)),
        foreign_entries=(
            ("sp|F1|FOREIGN_ONE foreign", "MKWVTFISLLFLFSSAYSR"),
            ("sp|F2|FOREIGN_TWO foreign", "AGGDDKYTFDKQVNTLPWR"),
        ),
    )

    assert result.n_entrapment == 1
    assert result.entrapment_pairs_path is None
    assert not tuple(tmp_path.glob("*.entrapment-pairs.parquet"))
