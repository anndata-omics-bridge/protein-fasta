"""Tests for generic FASTA header splitting."""

from protein_fasta.reading.header import ParsedHeader, parse_header


def test_parse_header_uses_first_nonempty_token_and_normalizes_description() -> None:
    assert parse_header(">  sp|P1|ONE   Protein kinase,   alpha ") == ParsedHeader(
        id="sp|P1|ONE",
        description="Protein kinase, alpha",
    )
    assert parse_header("P1") == ParsedHeader("P1", None)
    assert parse_header("   ") == ParsedHeader("", None)
