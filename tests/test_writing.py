"""Tests for generic FASTA serialization."""

from pathlib import Path

from protein_fasta.reading.parser import FastaRecord, read_records
from protein_fasta.reading.writer import write_records


def test_write_records_wraps_and_round_trips_lexical_values(tmp_path: Path) -> None:
    path = tmp_path / "written.fasta"
    records = [
        FastaRecord("P1 exact  header", "ABCDEFGHI"),
        FastaRecord("P2", ""),
    ]

    write_records(records, path, line_width=4)

    assert path.read_text() == ">P1 exact  header\nABCD\nEFGH\nI\n>P2\n"
    assert list(read_records(path)) == records


def test_write_records_can_disable_wrapping(tmp_path: Path) -> None:
    path = tmp_path / "unwrapped.fasta"
    write_records([FastaRecord("P1", "ABCDEF")], path, line_width=0)
    assert path.read_text() == ">P1\nABCDEF\n"
