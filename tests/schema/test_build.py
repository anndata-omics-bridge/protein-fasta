"""Tests for database-build storage documents."""

import datetime
from pathlib import Path

import pytest
from pydantic import ValidationError

from protein_fasta.database_build import (
    BuildDecoyOverride,
    DatabaseBuildOverrides,
    resolve_database_build,
)
from protein_fasta.schema.build import (
    DatabaseBuildDocument,
    DatabaseBuildProfileDocument,
    DatabaseBuildRequestDocument,
    DecoyDocument,
    DecoyMode,
    MetadataDocument,
    NamingDocument,
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


def _request(**updates: object) -> DatabaseBuildRequestDocument:
    values: dict[str, object] = {
        "targets": (Path("target.fasta"),),
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
        default_decoy=DecoyDocument(mode=DecoyMode.SHUFFLE, seed=11),
    )
    request = _request(
        metadata=MetadataDocument(org="request"),
        decoy=DecoyDocument(mode=DecoyMode.DECOYPYRAT, seed=17),
    )

    effective = resolve_database_build(
        profile,
        request,
        profile_base=profile_base,
        request_base=request_base,
        overrides=DatabaseBuildOverrides(
            date=datetime.date(2026, 8, 28),
            decoy=BuildDecoyOverride.REVERSE,
        ),
    )

    assert effective.metadata.org == "request"
    assert effective.date == datetime.date(2026, 8, 28)
    assert effective.decoy is not None
    assert effective.decoy.mode is DecoyMode.REVERSE
    assert effective.decoy.seed == 17
    assert effective.targets == ((request_base / "target.fasta").resolve(),)
    assert effective.output_dir == (request_base / "out").resolve()
    assert effective.diagnostics == (profile_base / "diagnostics.json").resolve()


def test_explicit_request_null_disables_profile_decoys(tmp_path: Path) -> None:
    profile = DatabaseBuildProfileDocument(
        default_decoy=DecoyDocument(mode=DecoyMode.SHUFFLE, seed=11)
    )
    request = _request(decoy=None)

    effective = resolve_database_build(
        profile,
        request,
        profile_base=tmp_path,
        request_base=tmp_path,
    )

    assert effective.decoy is None


def test_omitted_request_decoy_inherits_profile(tmp_path: Path) -> None:
    profile = DatabaseBuildProfileDocument(
        default_decoy=DecoyDocument(mode=DecoyMode.SHUFFLE, seed=11)
    )

    effective = resolve_database_build(
        profile,
        _request(),
        profile_base=tmp_path,
        request_base=tmp_path,
    )

    assert effective.decoy == profile.default_decoy


def test_cli_none_override_disables_request_decoys(tmp_path: Path) -> None:
    effective = resolve_database_build(
        DatabaseBuildProfileDocument(),
        _request(decoy=DecoyDocument(mode=DecoyMode.SHUFFLE)),
        profile_base=tmp_path,
        request_base=tmp_path,
        overrides=DatabaseBuildOverrides(decoy=BuildDecoyOverride.NONE),
    )

    assert effective.decoy is None
