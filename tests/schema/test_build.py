"""Tests for database-build storage documents."""

import datetime
from pathlib import Path

import pytest
from pydantic import ValidationError

from protein_fasta.schema.build import DatabaseBuildDocument, NamingDocument


def test_naming_rejects_unknown_template_variable() -> None:
    with pytest.raises(ValidationError, match="unsupported fields"):
        NamingDocument(dbname={"project": "p{project}_{surprise}"})


def test_naming_rejects_unsafe_template_expression() -> None:
    with pytest.raises(ValidationError, match="unsafe fields"):
        NamingDocument(dbname={"project": "p{project.__class__}"})


def test_naming_requires_every_filename_product() -> None:
    with pytest.raises(ValidationError, match="filename templates are missing"):
        NamingDocument(filename={"nondecoy": "{dbname}_{date}.{extension}"})


def test_build_rejects_unknown_selected_template(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="is not configured"):
        DatabaseBuildDocument(
            targets=(tmp_path / "target.fasta",),
            output_dir=tmp_path,
            date=datetime.date(2026, 8, 27),
            name_fields={"description": "demo"},
            template="unknown",
        )


def test_build_rejects_unknown_name_field(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="name_fields contains unsupported"):
        DatabaseBuildDocument(
            targets=(tmp_path / "target.fasta",),
            output_dir=tmp_path,
            date=datetime.date(2026, 8, 27),
            name_fields={"surprise": "demo"},
        )
