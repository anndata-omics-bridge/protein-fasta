"""Tests for biological-build storage documents and resolution."""

from __future__ import annotations

import datetime
import inspect
from pathlib import Path

import pytest
from pydantic import ValidationError

from protein_fasta.database_build import DatabaseBuildOverrides, resolve_database_build
from protein_fasta.schema.build import (
    DatabaseBuildProfileDocument,
    DatabaseBuildRequestDocument,
    EffectiveDatabaseBuildDocument,
    ForeignSpeciesEntrapmentDocument,
    MetadataDocument,
    NamingDocument,
    ShuffledEntrapmentDocument,
)


def test_naming_rejects_unknown_template_variable() -> None:
    with pytest.raises(ValidationError, match="unsupported fields"):
        NamingDocument(dbname={"project": "p{project}_{surprise}"})


def test_naming_rejects_unsafe_template_expression() -> None:
    with pytest.raises(ValidationError, match="unsafe fields"):
        NamingDocument(dbname={"project": "p{project.__class__}"})


def test_naming_requires_every_filename_product() -> None:
    with pytest.raises(ValidationError, match="filename templates are missing"):
        NamingDocument(filename={"nondecoy": "{dbname}_{date}.{extension}"})


def _request(**updates: object) -> DatabaseBuildRequestDocument:
    values: dict[str, object] = {
        "output_dir": Path("out"),
        "date": datetime.date(2026, 8, 27),
        "name_fields": {"description": "demo"},
    }
    values.update(updates)
    return DatabaseBuildRequestDocument.model_validate(values)


def test_build_resolution_records_profile_request_and_cli_precedence(tmp_path: Path) -> None:
    profile_base = tmp_path / "profiles"
    request_base = tmp_path / "requests"
    profile = DatabaseBuildProfileDocument(
        metadata=MetadataDocument(org="profile"),
        diagnostics=Path("diagnostics.json"),
    )
    request = _request(metadata=MetadataDocument(org="request"))

    effective = resolve_database_build(
        profile,
        request,
        profile_base=profile_base,
        request_base=request_base,
        overrides=DatabaseBuildOverrides(date=datetime.date(2026, 8, 28)),
    )

    assert effective.metadata.org == "request"
    assert effective.date == datetime.date(2026, 8, 28)
    assert effective.output_dir == (request_base / "out").resolve()
    assert effective.diagnostics == (profile_base / "diagnostics.json").resolve()


def test_effective_build_rejects_unknown_selected_template(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="is not configured"):
        EffectiveDatabaseBuildDocument(
            output_dir=tmp_path,
            date=datetime.date(2026, 8, 27),
            name_fields={"description": "demo"},
            template="unknown",
            naming=NamingDocument(),
            metadata=MetadataDocument(),
        )


def test_effective_build_rejects_unknown_name_field(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="name_fields contains unsupported"):
        EffectiveDatabaseBuildDocument(
            output_dir=tmp_path,
            date=datetime.date(2026, 8, 27),
            name_fields={"surprise": "demo"},
            template="project",
            naming=NamingDocument(),
            metadata=MetadataDocument(),
        )


def test_biological_build_documents_expose_no_decoy_policy() -> None:
    for model in (
        DatabaseBuildProfileDocument,
        DatabaseBuildRequestDocument,
        EffectiveDatabaseBuildDocument,
    ):
        assert "decoy" not in model.model_fields
    assert "decoy" not in inspect.signature(DatabaseBuildOverrides).parameters


def test_entrapment_variants_reject_fields_owned_by_the_other_strategy() -> None:
    with pytest.raises(ValidationError, match="reject_shared_foreign"):
        ShuffledEntrapmentDocument.model_validate(
            {"type": "shuffled", "reject_shared_foreign": False}
        )
    with pytest.raises(ValidationError, match="fix_peptide_n_term"):
        ForeignSpeciesEntrapmentDocument.model_validate(
            {"type": "foreign_species", "fix_peptide_n_term": False}
        )
