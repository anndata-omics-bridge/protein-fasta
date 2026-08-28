from __future__ import annotations

import datetime

from protein_fasta.build.metadata import build_sentinel_header
from protein_fasta.registry.metadata import parse_database_metadata
from protein_fasta.schema.build import MetadataDocument

CFG = MetadataDocument()


def test_parse_iso_dated_sentinel() -> None:
    header = (
        "aa|p42261_db1_9606wIsoforms_plus_custom|2026-06-05, customDB for p42261 complete "
        "uniprotkb of human including isoforms plus custom predicted ORFs, generated w prozor "
        "and installed by jonas, fgcz"
    )
    info = parse_database_metadata(header, CFG)
    assert info is not None
    assert info.dbname == "p42261_db1_9606wIsoforms_plus_custom"
    assert info.date == "2026-06-05"
    assert info.description.startswith("customDB for p42261")


def test_parse_compact_dated_sentinel() -> None:
    header = (
        "aa|p36602_db1_SwissProt_TrEMBL_Mleprae|20251210 SwissProt (20250403) AND TrEMBL "
        "for Mycobacterium leprae (ID:1769, downloaded 20250818)"
    )
    info = parse_database_metadata(header, CFG)
    assert info is not None
    assert info.dbname == "p36602_db1_SwissProt_TrEMBL_Mleprae"
    assert info.date == "20251210"
    assert info.description.startswith("SwissProt")


def test_section_marker_is_not_a_db_sentinel() -> None:
    header = "aa|Cont_UniversalContaminants|fgcz_universal_contaminants_github_20241112 file downloaded from ..."
    assert parse_database_metadata(header, CFG) is None


def test_non_sentinel_returns_none() -> None:
    assert parse_database_metadata("sp|P12345|FOO_HUMAN OS=Homo sapiens", CFG) is None


def test_leading_gt_is_tolerated() -> None:
    info = parse_database_metadata(">aa|p1_db1|2026-07-01 desc", CFG)
    assert info is not None
    assert info.dbname == "p1_db1"


def test_build_parse_roundtrip() -> None:
    header = build_sentinel_header(
        "p999_db1_test", "my description", datetime.date(2026, 7, 1), CFG
    )
    info = parse_database_metadata(header, CFG)
    assert info is not None
    assert info.dbname == "p999_db1_test"
    assert info.date == "2026-07-01"
    assert info.description.startswith("my description")
