"""Tests for strict JSON document loading."""

import json
from importlib import resources
from pathlib import Path
from typing import cast

import pytest
from pydantic import BaseModel, ValidationError

from protein_fasta.documents import (
    load_builtin_diagnostic_document,
    load_builtin_entry_classifier_document,
    load_diagnostic_document,
)
from protein_fasta.schema.analytics import DigestionDocument, EnzymeDocument
from protein_fasta.schema.build import DatabaseBuildDocument
from protein_fasta.schema.diagnostics import (
    DiagnosticDocument,
    EntryClassifierCatalogDocument,
    EntryClassifierDocument,
)
from protein_fasta.schema.frame_formats import HeaderFormatDocument
from protein_fasta.schema.registry import RegistryDiagnosticDocument, RegistryDocument


def test_builtin_documents_load_with_stable_versions() -> None:
    assert load_builtin_diagnostic_document().file_version == "1"
    assert load_builtin_entry_classifier_document().file_version == "2"


def test_explicit_loader_names_invalid_source_path(tmp_path: Path) -> None:
    path = tmp_path / "diagnostics.json"
    path.write_text('{"schema_version": "0.1"}')

    with pytest.raises(ValueError, match=r"diagnostics\.json"):
        load_diagnostic_document(path)


@pytest.mark.parametrize(
    ("field", "pattern"),
    [
        ("removable_prefix_patterns", "REV_"),
        ("removable_prefix_patterns", "^"),
        ("removable_suffix_patterns", "_target"),
        ("removable_suffix_patterns", "$"),
    ],
)
def test_classifier_rejects_unsafe_removable_patterns(field: str, pattern: str) -> None:
    values: dict[str, object] = {
        "name": "bad",
        "output_column": "is_bad",
        field: [pattern],
    }
    with pytest.raises(ValidationError):
        EntryClassifierDocument.model_validate(values)


def test_classifier_requires_at_least_one_pattern() -> None:
    with pytest.raises(ValidationError, match="at least one pattern"):
        EntryClassifierDocument(name="empty", output_column="is_empty")


@pytest.mark.parametrize(
    ("name", "model"),
    [
        ("database_build.schema.json", DatabaseBuildDocument),
        ("diagnostic.schema.json", DiagnosticDocument),
        ("digestion.schema.json", DigestionDocument),
        ("entry_classifier.schema.json", EntryClassifierCatalogDocument),
        ("enzyme.schema.json", EnzymeDocument),
        ("header_format.schema.json", HeaderFormatDocument),
        ("registry.schema.json", RegistryDocument),
        ("registry_diagnostic.schema.json", RegistryDiagnosticDocument),
    ],
)
def test_committed_json_schema_matches_pydantic(
    name: str,
    model: type[BaseModel],
) -> None:
    resource = resources.files("protein_fasta").joinpath("documents", "_schema", name)
    committed = cast(
        dict[str, object],
        json.loads(resource.read_text(encoding="utf-8")),
    )
    assert committed == model.model_json_schema()
