"""Tests for streaming FASTA input."""

from __future__ import annotations

import bz2
import gzip
from io import StringIO
from pathlib import Path

import pytest

from protein_fasta.reading.parser import parse_records, parse_text, read_headers, read_records
from protein_fasta.reading.record import FastaReadError, FastaRecord


@pytest.mark.parametrize("suffix", [".fasta", ".fasta.gz", ".fasta.bz2"])
def test_read_records_streams_plain_and_compressed_sources(tmp_path: Path, suffix: str) -> None:
    path = tmp_path / f"proteins{suffix}"
    text = ">sp|P1|ONE first\nAA C\nD*\n\n>sp|P2|TWO\nEF-G\n"
    if suffix.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    elif suffix.endswith(".bz2"):
        with bz2.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text)

    assert list(read_records(path)) == [
        FastaRecord("sp|P1|ONE first", "AACD*"),
        FastaRecord("sp|P2|TWO", "EF-G"),
    ]
    assert list(read_headers(path)) == ["sp|P1|ONE first", "sp|P2|TWO"]


def test_parse_text_is_explicit_and_preserves_empty_records() -> None:
    assert list(parse_text(">empty\n>full\nPEP\n")) == [
        FastaRecord("empty", ""),
        FastaRecord("full", "PEP"),
    ]


def test_parse_records_rejects_content_before_header_with_location() -> None:
    with pytest.raises(FastaReadError) as caught:
        list(parse_records(StringIO("SEQUENCE\n"), source_name="upload"))

    assert caught.value.source_name == "upload"
    assert caught.value.line_number == 1
    assert str(caught.value) == ("upload: sequence content before the first FASTA header at line 1")


def test_parse_records_accepts_empty_and_blank_input() -> None:
    assert list(parse_records([])) == []
    assert list(parse_records(["\n", " \t\n"])) == []


def test_read_records_reports_non_utf8_path(tmp_path: Path) -> None:
    path = tmp_path / "legacy.fasta"
    path.write_bytes(b">P1\n\xff\n")

    with pytest.raises(FastaReadError, match=r"legacy\.fasta: input is not valid UTF-8"):
        list(read_records(path))


def test_read_records_reports_missing_path(tmp_path: Path) -> None:
    path = tmp_path / "missing.fasta"

    with pytest.raises(FastaReadError, match=r"missing\.fasta: input cannot be read"):
        list(read_records(path))

    with pytest.raises(FastaReadError, match=r"missing\.fasta: input cannot be read"):
        list(read_headers(path))
