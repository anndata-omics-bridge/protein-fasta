"""Tests for strict frame-format storage documents and compilation."""

import pytest
from pydantic import ValidationError

from protein_fasta.frame_compile import make_frame_parsers
from protein_fasta.schema.frame_formats import (
    HeaderColumnDocument,
    HeaderColumnsDocument,
    HeaderFormatCatalogDocument,
    HeaderFormatDocument,
)


@pytest.mark.parametrize(
    "values",
    [
        {"name": "bad", "type": "string"},
        {"name": "bad", "pattern": "(.+)", "value": "both"},
        {"name": "bad", "pattern": "no capture"},
        {"name": "bad", "pattern": "(.+)(extra)"},
        {"name": "id", "value": "replacement"},
        {"name": "bad", "type": "integer", "pattern": "(.+)", "values": {"1": "one"}},
    ],
)
def test_column_document_rejects_invalid_sources(values: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        HeaderColumnDocument.model_validate(values)


def test_polars_compilation_rejects_unsupported_lookahead() -> None:
    document = HeaderFormatDocument(
        file_version="1",
        format="lookahead",
        detection_pattern=r"^P(?=1)",
        columns=HeaderColumnsDocument(),
    )

    with pytest.raises(ValueError, match="not supported by Polars"):
        make_frame_parsers(HeaderFormatCatalogDocument(formats=(document,)))
