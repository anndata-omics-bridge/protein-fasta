from __future__ import annotations

from dataclasses import fields

import polars as pl
import pytest

from protein_fasta.analytics.peptide_properties import PeptideProperties
from protein_fasta.peptide_frame import peptide_property_frame


def test_one_row_per_sequence_in_input_order() -> None:
    frame = peptide_property_frame(["KQPWWR", "EAAAMGPTK", "KQPWWR"])

    assert frame.get_column("sequence").to_list() == ["KQPWWR", "EAAAMGPTK", "KQPWWR"]
    assert frame.get_column("length").to_list() == [6, 9, 6]
    assert frame.get_column("c_terminal_residue").to_list() == ["R", "K", "R"]


def test_schema_follows_the_property_fields() -> None:
    frame = peptide_property_frame([])

    assert frame.height == 0
    assert frame.columns == ["sequence", *(field.name for field in fields(PeptideProperties))]
    assert frame.schema["length"] == pl.Int64
    assert frame.schema["isoelectric_point"] == pl.Float64
    assert frame.schema["c_terminal_residue"] == pl.String
    assert frame.schema["contains_ng_motif"] == pl.Boolean


def test_non_standard_sequence_has_null_properties() -> None:
    frame = peptide_property_frame(["RSGASILQAGCUG", "KQPWWR"])

    assert frame.get_column("sequence").to_list() == ["RSGASILQAGCUG", "KQPWWR"]
    assert frame.drop("sequence").row(0) == (None,) * (frame.width - 1)
    assert frame.get_column("length").to_list() == [None, 6]


def test_malformed_sequence_is_rejected() -> None:
    with pytest.raises(ValueError, match="upper-case letters"):
        peptide_property_frame(["KQPWWR", ""])
