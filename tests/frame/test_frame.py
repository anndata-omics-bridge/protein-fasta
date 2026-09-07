"""Tests for the base and configured Polars frame APIs."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from protein_fasta.analytics.hashing import file_checksum
from protein_fasta.documents import (
    load_builtin_entry_classifier_document,
    load_builtin_header_format_catalog,
)
from protein_fasta.frame import (
    ProteinDatabase,
    ProteinFormat,
    read_basic_protein_frame,
    read_configured_protein_frame,
    read_header_format_diagnostics_frame,
    read_protein_frame,
    read_strict_configured_protein_frame,
    read_strict_protein_frame,
    refseq,
    uniprotkb,
)
from protein_fasta.frame_formats.extraction import FrameExtractionError
from protein_fasta.schema.diagnostics import (
    EntryClassifierCatalogDocument,
    EntryClassifierDocument,
)
from protein_fasta.schema.frame_formats import HeaderFormatCatalogDocument


def _write_fasta(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "proteins.fasta"
    path.write_text(text)
    return path


def test_basic_frame_has_exact_stable_schema_and_normalized_values(tmp_path: Path) -> None:
    path = _write_fasta(tmp_path, ">P1   description  here\nac d*\n>P2\nEF\n")

    frame = read_basic_protein_frame(path)

    assert frame.schema == {
        "id": pl.String,
        "description": pl.String,
        "sequence": pl.String,
    }
    assert frame.to_dicts() == [
        {"id": "P1", "description": "description here", "sequence": "ACD"},
        {"id": "P2", "description": None, "sequence": "EF"},
    ]


def test_empty_frame_has_exact_base_schema(tmp_path: Path) -> None:
    path = _write_fasta(tmp_path, "")
    assert read_protein_frame(path).schema == {
        "id": pl.String,
        "description": pl.String,
        "sequence": pl.String,
    }


def test_uniprot_frame_peels_decorations_and_extracts_best_columns(tmp_path: Path) -> None:
    path = _write_fasta(
        tmp_path,
        (
            ">REV_CON__sp|P27748|ACOX_RALEH Acetoin catabolism protein X "
            "OS=Ralstonia eutropha OX=10116 GN=acoX PE=4 SV=2\nacd*\n"
            ">tr|Q00001|RL40_YEAST Ribosomal protein OS=Saccharomyces cerevisiae "
            "OX=559292 PE=3 SV=1\nEFG\n"
        ),
    )

    frame = read_protein_frame(path)

    assert frame.columns == [
        "id",
        "description",
        "sequence",
        "is_decoy",
        "is_contaminant",
        "database",
        "review_status",
        "accession",
        "entry_name",
        "entry_mnemonic",
        "organism_mnemonic",
        "protein_name",
        "organism_name",
        "taxonomy_id",
        "gene_name",
        "protein_existence",
        "sequence_version",
    ]
    first, second = frame.to_dicts()
    assert first["id"] == "REV_CON__sp|P27748|ACOX_RALEH"
    assert first["is_decoy"] is True
    assert first["is_contaminant"] is True
    assert first["database"] == "uniprotkb"
    assert first["review_status"] == "reviewed"
    assert first["entry_name"] == "ACOX_RALEH"
    assert first["entry_mnemonic"] == "ACOX"
    assert first["organism_mnemonic"] == "RALEH"
    assert first["organism_name"] == "Ralstonia eutropha"
    assert first["taxonomy_id"] == 10116
    assert first["gene_name"] == "acoX"
    assert second["review_status"] == "unreviewed"
    assert second["gene_name"] is None


def test_configured_suffix_is_removed_before_uniprot_parsing(tmp_path: Path) -> None:
    path = _write_fasta(
        tmp_path,
        ">sp|P12345|RL40_YEAST_p_target Protein OS=Yeast OX=559292 PE=1 SV=1\nAA\n",
    )

    classifiers = EntryClassifierCatalogDocument(
        file_version="test",
        classifiers=(
            EntryClassifierDocument(
                name="entrapment",
                output_column="is_entrapment",
                removable_suffix_patterns=("_p_target$",),
            ),
        ),
    )
    row = read_configured_protein_frame(
        path,
        load_builtin_header_format_catalog(),
        classifiers,
    ).row(0, named=True)

    assert row["id"] == "sp|P12345|RL40_YEAST_p_target"
    assert row["is_entrapment"] is True
    assert row["entry_name"] == "RL40_YEAST"
    assert row["organism_mnemonic"] == "YEAST"


def test_refseq_frame_extracts_accession_name_and_optional_organism(tmp_path: Path) -> None:
    path = _write_fasta(
        tmp_path,
        (
            ">ref|NP_123456.1| ATP synthase subunit [Homo sapiens]\nAA\n"
            ">ref|WP_123456789.1| Hypothetical protein\nBB\n"
        ),
    )

    frame = read_protein_frame(path)

    assert frame["database"].to_list() == ["refseq", "refseq"]
    assert frame["accession"].to_list() == ["NP_123456.1", "WP_123456789.1"]
    assert frame["protein_name"].to_list() == [
        "ATP synthase subunit",
        "Hypothetical protein",
    ]
    assert frame["organism_name"].to_list() == ["Homo sapiens", None]


def test_mixed_uniprot_and_refseq_rows_are_enriched_independently(
    tmp_path: Path,
) -> None:
    path = _write_fasta(
        tmp_path,
        (
            ">sp|P12345|RL40_YEAST Protein OS=Yeast OX=1 PE=1 SV=1\nAA\n"
            ">ref|NP_123456.1| ATP synthase subunit [Homo sapiens]\nBB\n"
        ),
    )
    classifiers = load_builtin_entry_classifier_document()
    diagnostics = read_header_format_diagnostics_frame(
        path,
        load_builtin_header_format_catalog(),
        classifiers,
    )

    assert diagnostics.to_dicts() == [
        {
            "format": "refseq",
            "matched_rows": 1,
            "total_rows": 2,
            "status": "partial",
        },
        {
            "format": "uniprotkb",
            "matched_rows": 1,
            "total_rows": 2,
            "status": "partial",
        },
    ]
    frame = read_configured_protein_frame(
        path,
        load_builtin_header_format_catalog(),
        classifiers,
    )
    uniprot, refseq = frame.to_dicts()
    assert uniprot["database"] == "uniprotkb"
    assert uniprot["accession"] == "P12345"
    assert uniprot["review_status"] == "reviewed"
    assert refseq["database"] == "refseq"
    assert refseq["accession"] == "NP_123456.1"
    assert refseq["protein_name"] == "ATP synthase subunit"
    assert refseq["review_status"] is None
    assert read_protein_frame(path).to_dicts() == frame.to_dicts()


def test_strict_mixed_frame_returns_exact_base_columns(tmp_path: Path) -> None:
    path = _write_fasta(
        tmp_path,
        (
            ">sp|P12345|RL40_YEAST Protein OS=Yeast OX=1 PE=1 SV=1\nAA\n"
            ">ref|NP_123456.1| ATP synthase subunit [Homo sapiens]\nBB\n"
        ),
    )
    classifiers = load_builtin_entry_classifier_document()

    assert read_strict_protein_frame(path).columns == ["id", "description", "sequence"]
    assert read_strict_configured_protein_frame(
        path,
        load_builtin_header_format_catalog(),
        classifiers,
    ).columns == ["id", "description", "sequence"]


def test_recognized_rows_are_enriched_beside_unknown_rows(tmp_path: Path) -> None:
    frame = read_protein_frame(
        _write_fasta(
            tmp_path,
            ">sp|P12345|RL40_YEAST Protein OS=Yeast OX=1 PE=1 SV=1\nAA\n>P1 generic\nBB\n",
        )
    )

    assert frame["database"].to_list() == ["uniprotkb", None]
    assert frame["id"].to_list() == ["sp|P12345|RL40_YEAST", "P1"]
    assert frame["sequence"].to_list() == ["AA", "BB"]


def test_wholly_unknown_frame_returns_exact_base_columns(tmp_path: Path) -> None:
    assert read_protein_frame(_write_fasta(tmp_path, ">P1 generic\nAA\n")).columns == [
        "id",
        "description",
        "sequence",
    ]


def test_diagnostics_frame_reports_ambiguous_complete_matches(tmp_path: Path) -> None:
    path = _write_fasta(tmp_path, ">sp|P12345|RL40_YEAST Protein\nAA\n")
    built_in = next(
        document
        for document in load_builtin_header_format_catalog().formats
        if document.format == "uniprotkb"
    )
    catalog = HeaderFormatCatalogDocument(
        formats=(built_in, built_in.model_copy(update={"format": "same_shape"}))
    )

    diagnostics = read_header_format_diagnostics_frame(
        path,
        catalog,
        load_builtin_entry_classifier_document(),
    )

    assert diagnostics["status"].to_list() == ["ambiguous", "ambiguous"]
    assert read_configured_protein_frame(
        path,
        catalog,
        load_builtin_entry_classifier_document(),
    ).columns == ["id", "description", "sequence"]


def test_selected_format_raises_when_required_extraction_is_missing(tmp_path: Path) -> None:
    path = _write_fasta(tmp_path, ">sp|P12345|RL40_YEAST Protein\nAA\n")
    built_in = next(
        document
        for document in load_builtin_header_format_catalog().formats
        if document.format == "uniprotkb"
    )
    required = built_in.columns.required[1].model_copy(update={"pattern": r"^missing=(.+)$"})
    broken = built_in.model_copy(
        update={"columns": built_in.columns.model_copy(update={"required": (required,)})}
    )

    with pytest.raises(FrameExtractionError, match="review_status"):
        read_configured_protein_frame(
            path,
            HeaderFormatCatalogDocument(formats=(broken,)),
            load_builtin_entry_classifier_document(),
        )


def test_protein_database_combines_sources_with_stable_provenance(tmp_path: Path) -> None:
    human = tmp_path / "human.fasta"
    human.write_text(
        ">sp|P12345|RL40_YEAST Protein OS=Yeast OX=1 PE=1 SV=1\nAA\n",
        encoding="utf-8",
    )
    contaminants = tmp_path / "contaminants.fasta"
    contaminants.write_text(
        ">ref|NP_123456.1| ATP synthase subunit [Homo sapiens]\nBB\n",
        encoding="utf-8",
    )

    proteins = ProteinDatabase(uniprotkb, refseq).parse((human, contaminants))

    assert proteins["id"].to_list() == [
        "sp|P12345|RL40_YEAST",
        "ref|NP_123456.1|",
    ]
    assert proteins["database"].to_list() == ["uniprotkb", "refseq"]
    assert proteins["fasta_source_path"].to_list() == [str(human), str(contaminants)]
    assert proteins["fasta_source_checksum"].to_list() == [
        file_checksum(human),
        file_checksum(contaminants),
    ]
    assert proteins["fasta_source_ordinal"].to_list() == [0, 1]
    assert proteins["fasta_record_ordinal"].to_list() == [0, 0]


def test_protein_database_empty_parse_has_configured_schema() -> None:
    proteins = ProteinDatabase(uniprotkb, refseq).parse(())

    assert proteins.is_empty()
    assert proteins.columns[:5] == [
        "id",
        "description",
        "sequence",
        "is_decoy",
        "is_contaminant",
    ]
    assert proteins.columns[-4:] == [
        "fasta_source_path",
        "fasta_source_checksum",
        "fasta_source_ordinal",
        "fasta_record_ordinal",
    ]


@pytest.mark.parametrize(
    ("formats", "message"),
    [
        ((), "at least one"),
        ((uniprotkb, uniprotkb), "must be unique"),
        ((ProteinFormat("missing"),), "unknown packaged protein format"),
    ],
)
def test_protein_database_rejects_invalid_formats(
    formats: tuple[ProteinFormat, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        ProteinDatabase(*formats)
