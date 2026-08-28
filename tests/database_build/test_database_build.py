"""Shared database assembly keeps construction separate from registry analysis."""

from __future__ import annotations

import datetime
from pathlib import Path

import pytest

from protein_fasta.database_build import (
    ConflictingIdError,
    ContaminantBlock,
    Entry,
    PipelineResult,
    build_database,
)
from protein_fasta.reading.header import parse_header
from protein_fasta.reading.parser import read_records
from protein_fasta.registry.rules import load_registry_diagnostics
from protein_fasta.schema.build import (
    EntrapmentDocument,
    EntrapmentStrategy,
    MetadataDocument,
    NamingDocument,
)


def _build(
    tmp_path: Path,
    *,
    targets: tuple[Entry, ...] = (("sp|P1|ONE Protein one", "ak*"),),
    contaminant_blocks: tuple[ContaminantBlock, ...] = (),
) -> PipelineResult:
    return build_database(
        targets=targets,
        name_fields={"description": "demo"},
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
        diagnostics=load_registry_diagnostics(),
        output_dir=tmp_path,
        date=datetime.date(2026, 8, 27),
        template="derived",
        contaminant_blocks=contaminant_blocks,
        add_decoys=False,
    )


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


def test_entrapment_joins_the_biological_space_before_decoys(tmp_path: Path) -> None:
    diagnostics = load_registry_diagnostics()
    result = build_database(
        targets=(
            ("sp|P1|ONE first", "MPEPTIDEKAAALEVSSKMWNQTR"),
            ("sp|P2|TWO second", "MSAMPLERGGGYTFDKQVNTLPWR"),
        ),
        name_fields={"description": "entrapment"},
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
        diagnostics=diagnostics,
        output_dir=tmp_path,
        date=datetime.date(2026, 8, 27),
        template="derived",
        add_decoys=True,
        entrapment_spec=EntrapmentDocument(fold=1, seed=7),
    )

    identifiers = [parse_header(record.raw_header).id for record in read_records(result.path)]
    entrapment_ids = {
        identifier
        for identifier in identifiers
        if identifier.endswith("_p_target") and not identifier.startswith(diagnostics.decoy_prefix)
    }

    assert result.n_entrapment == 2
    assert result.n_decoy == result.n_target + result.n_contaminant + result.n_entrapment
    assert len(entrapment_ids) == result.n_entrapment
    assert {f"REV_{identifier}" for identifier in entrapment_ids} <= set(identifiers)
    assert min(identifiers.index(identifier) for identifier in entrapment_ids) < min(
        identifiers.index(f"REV_{identifier}") for identifier in entrapment_ids
    )
    assert result.entrapment_pairs_path is not None
    assert result.entrapment_pairs_path.is_file()
    assert result.entrapment_strategy == EntrapmentStrategy.SHUFFLED
    assert result.entrapment_requested_fold == result.entrapment_achieved_fold == 1


def test_foreign_species_entrapment_does_not_claim_peptide_pairs(tmp_path: Path) -> None:
    result = build_database(
        targets=(("sp|P1|ONE first", "MPEPTIDEKTESTSEQUENCEK"),),
        name_fields={"description": "foreign"},
        naming=NamingDocument(default_dbname="derived"),
        metadata=MetadataDocument(),
        diagnostics=load_registry_diagnostics(),
        output_dir=tmp_path,
        date=datetime.date(2026, 8, 27),
        template="derived",
        add_decoys=False,
        entrapment_spec=EntrapmentDocument(
            strategy=EntrapmentStrategy.FOREIGN_SPECIES,
            fold=1,
        ),
        foreign_entries=(
            ("sp|F1|FOREIGN_ONE foreign", "MKWVTFISLLFLFSSAYSR"),
            ("sp|F2|FOREIGN_TWO foreign", "AGGDDKYTFDKQVNTLPWR"),
        ),
    )

    assert result.n_entrapment == 1
    assert result.entrapment_pairs_path is None
    assert not tuple(tmp_path.glob("*.entrapment_pairs.tsv"))
