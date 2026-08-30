from __future__ import annotations

import datetime

import pytest

from protein_fasta.database.naming import build_dbname, build_description, build_fasta_name
from protein_fasta.database_compile import make_database_naming
from protein_fasta.registry.filenames import parse_filename
from protein_fasta.schema.build import NamingDocument

CFG_DOCUMENT = NamingDocument()
CFG = make_database_naming(CFG_DOCUMENT)


def test_build_basic() -> None:
    name = build_fasta_name(
        config=CFG,
        template="project",
        project=41977,
        dbn=1,
        description="canisLupus",
        date=datetime.date(2026, 6, 3),
        decoy=False,
    )
    assert name == "p41977_db1_canisLupus_20260603.fasta"


def test_build_decoy() -> None:
    name = build_fasta_name(
        config=CFG,
        template="project",
        project=39154,
        dbn=1,
        description="UP000001449_Thalassiosira",
        date=datetime.date(2025, 8, 14),
        decoy=True,
    )
    assert name == "p39154_db1_UP000001449_Thalassiosira_d_20250814.fasta"


def test_default_dbname_is_project() -> None:
    assert build_dbname(config=CFG, project=1, dbn=2, description="x") == "p1_db2_x"


def test_fgcz_dbname() -> None:
    assert (
        build_dbname(config=CFG, template="fgcz", taxid=9606, description="reviewed")
        == "fgcz_9606_reviewed"
    )


def test_empty_field_collapses_separators() -> None:
    # A missing description must not leave a trailing/double separator.
    assert (
        build_dbname(config=CFG, template="project", project=999, dbn=1, description="")
        == "p999_db1"
    )


def test_unknown_template_raises() -> None:
    with pytest.raises(KeyError):
        build_dbname(config=CFG, template="nope", project=1, dbn=1, description="x")


@pytest.mark.parametrize(
    ("filename", "dbname", "decoy", "date"),
    [
        (
            "p41977_db1_canisLupus_20260603.fasta",
            "p41977_db1_canisLupus",
            False,
            datetime.date(2026, 6, 3),
        ),
        (
            "p39154_db1_UP000001449_Thalassiosira_d_20250814.fasta",
            "p39154_db1_UP000001449_Thalassiosira",
            True,
            datetime.date(2025, 8, 14),
        ),
        ("R_fgcz_6239_20050304.fasta", "fgcz_6239", True, datetime.date(2005, 3, 4)),
        (
            "fgcz_6239_wormpep200_20090812_d.fasta",
            "fgcz_6239_wormpep200",
            True,
            datetime.date(2009, 8, 12),
        ),
        # glued-date source dump: no standalone 8-digit token, so date is not recovered
        ("uniprot_sprot20240124.fasta", "uniprot_sprot20240124", False, None),
    ],
)
def test_parse_filename(
    filename: str, dbname: str, decoy: bool, date: datetime.date | None
) -> None:
    parsed = parse_filename(filename, CFG_DOCUMENT)
    assert parsed.dbname == dbname
    assert parsed.is_decoy is decoy
    assert parsed.date == date


def test_filename_patterns_are_config_driven() -> None:
    # Editing the filename pattern in config changes the output with no code change.
    cfg_document = NamingDocument(
        filename={
            **CFG.filename,
            "decoy": "{dbname}.decoy.{date}.{extension}",
            "nondecoy": "{dbname}.{date}.{extension}",
        }
    )
    cfg = make_database_naming(cfg_document)
    name = build_fasta_name(
        config=cfg,
        template="project",
        project=1,
        dbn=1,
        description="x",
        date=datetime.date(2026, 1, 2),
        decoy=True,
    )
    assert name == "p1_db1_x.decoy.20260102.fasta"


def test_build_description_joins_taxids_and_mode_flag() -> None:
    # A merge of Human + Yeast Swiss-Prot -> terse "9606_4932_sp".
    assert build_description(config=CFG, taxids=[9606, 4932], mode="swissprot") == "9606_4932_sp"


def test_build_description_empty_flag_collapses() -> None:
    # swissprot_trembl has an empty mode flag -> taxids only.
    assert build_description(config=CFG, taxids=[9606], mode="swissprot_trembl") == "9606"


def test_build_description_flag_only_when_no_taxids() -> None:
    # When the template already carries the taxid, the app passes no taxids.
    assert build_description(config=CFG, taxids=[], mode="one_seq_per_gene") == "1spg"


def test_build_parse_roundtrip() -> None:
    date = datetime.date(2026, 7, 1)
    name = build_fasta_name(
        config=CFG,
        template="project",
        project=999,
        dbn=2,
        description="human_reviewed",
        date=date,
        decoy=True,
    )
    parsed = parse_filename(name, CFG_DOCUMENT)
    assert parsed.dbname == "p999_db2_human_reviewed"
    assert parsed.is_decoy is True
    assert parsed.date == date


def test_build_parse_roundtrip_honours_configured_extension() -> None:
    # parse_filename must strip the configured extension, not just the legacy set.
    cfg_document = NamingDocument(extension="pep")
    cfg = make_database_naming(cfg_document)
    date = datetime.date(2026, 1, 2)
    name = build_fasta_name(
        config=cfg, template="project", project=1, dbn=1, description="x", date=date, decoy=False
    )
    assert name == "p1_db1_x_20260102.pep"
    parsed = parse_filename(name, cfg_document)
    assert parsed.dbname == "p1_db1_x"
    assert parsed.date == date
