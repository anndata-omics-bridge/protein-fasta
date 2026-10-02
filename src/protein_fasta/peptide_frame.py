"""Polars peptide-property frame API."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import astuple, fields

import polars as pl

from protein_fasta.analytics.peptide_properties import PeptideProperties, peptide_properties

# Field annotations are postponed, so `field.type` is the annotation text.
_POLARS_TYPES: dict[str, type[pl.DataType]] = {
    "int": pl.Int64,
    "float": pl.Float64,
    "str": pl.String,
    "bool": pl.Boolean,
}
_SCHEMA: dict[str, type[pl.DataType]] = {
    "sequence": pl.String,
    **{field.name: _POLARS_TYPES[str(field.type)] for field in fields(PeptideProperties)},
}
_NULL_PROPERTIES = (None,) * len(fields(PeptideProperties))


def peptide_property_frame(sequences: Iterable[str], /) -> pl.DataFrame:
    """Return one row of sequence-derived properties per peptide, in input order.

    Args:
        sequences: Stripped, upper-case peptide sequences; duplicates are kept.

    Returns:
        A frame with ``sequence`` followed by one column per ``PeptideProperties`` field. A
        sequence containing a residue other than the 20 standard amino acids has null properties.

    Raises:
        ValueError: If a sequence is empty or contains anything but upper-case ASCII letters.
    """
    rows: list[tuple[object, ...]] = []
    for sequence in sequences:
        properties = peptide_properties(sequence)
        rows.append((sequence, *(_NULL_PROPERTIES if properties is None else astuple(properties))))
    return pl.DataFrame(rows, schema=_SCHEMA, orient="row")
