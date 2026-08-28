"""Vectorized entry classification and identifier decoration peeling."""

from __future__ import annotations

import polars as pl

from protein_fasta.frame_formats.runtime import (
    CompiledFrameClassifiers,
    concat_strings,
    with_columns,
)

_ID = "id"
_DESCRIPTION = "description"
_WORKING_IDENTIFIER = "__working_identifier"
WORKING_HEADER = "__working_header"


def append_classifications(
    frame: pl.DataFrame,
    classifiers: CompiledFrameClassifiers,
    /,
) -> pl.DataFrame:
    """Append independent labels and the reconstructed working header."""
    working = frame.get_column(_ID)
    flags = {
        classifier.name: pl.Series(
            classifier.output_column,
            [False] * frame.height,
            dtype=pl.Boolean,
        )
        for classifier in classifiers.classifiers
    }

    while working.len() > 0:
        before_pass = working
        for classifier in classifiers.classifiers:
            for pattern in classifier.removable_prefix_patterns:
                matches = working.str.contains(pattern)
                flags[classifier.name] = flags[classifier.name] | matches
                working = working.str.replace(pattern, "", n=1)
            for pattern in classifier.removable_suffix_patterns:
                matches = working.str.contains(pattern)
                flags[classifier.name] = flags[classifier.name] | matches
                working = working.str.replace(pattern, "", n=1)
        if working.equals(before_pass):
            break

    for classifier in classifiers.classifiers:
        for pattern in classifier.match_patterns:
            flags[classifier.name] = flags[classifier.name] | working.str.contains(pattern)

    expressions: list[pl.Series] = [working.rename(_WORKING_IDENTIFIER)]
    expressions.extend(
        flags[classifier.name].rename(classifier.output_column)
        for classifier in classifiers.classifiers
    )
    classified = with_columns(frame, expressions)
    working_header = concat_strings(
        [pl.col(_WORKING_IDENTIFIER), pl.col(_DESCRIPTION)],
        separator=" ",
        ignore_nulls=True,
    ).alias(WORKING_HEADER)
    return with_columns(classified, [working_header]).drop(_WORKING_IDENTIFIER)
