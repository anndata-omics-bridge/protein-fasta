"""Load validated protein-FASTA JSON documents from paths and package data."""

from __future__ import annotations

import json
from importlib import resources
from pathlib import Path

from pydantic import BaseModel, ValidationError

from protein_fasta.schema.analytics import DigestionDocument, EnzymeDocument
from protein_fasta.schema.build import (
    DatabaseBuildProfileDocument,
    DatabaseBuildRequestDocument,
)
from protein_fasta.schema.candidate import CandidateRequestDocument
from protein_fasta.schema.decoy import DecoyRequestDocument
from protein_fasta.schema.decoy_report import DecoyReportRequestDocument
from protein_fasta.schema.diagnostics import (
    DiagnosticDocument,
    EntryClassifierCatalogDocument,
)
from protein_fasta.schema.frame_formats import (
    HeaderFormatCatalogDocument,
    HeaderFormatDocument,
)
from protein_fasta.schema.peptide import (
    PeptideBuildRequestDocument,
    PeptideComparisonRequestDocument,
)
from protein_fasta.schema.protein_input import (
    DerivedProteinInputRequestDocument,
    ProteinInputRequestDocument,
)
from protein_fasta.schema.registry import RegistryDocument
from protein_fasta.schema.uniprot import (
    UniProtCatalogRequestDocument,
    UniProtDownloadRequestDocument,
)

_PACKAGE_DOCUMENTS = resources.files("protein_fasta").joinpath("documents")


def load_diagnostic_document(path: Path, /) -> DiagnosticDocument:
    """Load one diagnostic document from an explicit JSON path."""
    return _load_path(path, DiagnosticDocument)


def load_entry_classifier_document(path: Path, /) -> EntryClassifierCatalogDocument:
    """Load one entry-classifier document from an explicit JSON path."""
    return _load_path(path, EntryClassifierCatalogDocument)


def load_builtin_diagnostic_document() -> DiagnosticDocument:
    """Load the packaged application-neutral diagnostic rules."""
    resource = _PACKAGE_DOCUMENTS.joinpath("diagnostics", "rules.json")
    return _load_text(resource.read_text(encoding="utf-8"), str(resource), DiagnosticDocument)


def load_builtin_entry_classifier_document() -> EntryClassifierCatalogDocument:
    """Load the packaged common entry decorations."""
    resource = _PACKAGE_DOCUMENTS.joinpath("entry_classifiers", "rules.json")
    return _load_text(
        resource.read_text(encoding="utf-8"),
        str(resource),
        EntryClassifierCatalogDocument,
    )


def load_enzyme_document(path: Path, /) -> EnzymeDocument:
    """Load one enzyme rule from an explicit JSON path."""
    return _load_path(path, EnzymeDocument)


def load_digestion_document(path: Path, /) -> DigestionDocument:
    """Load one digestion configuration from an explicit JSON path."""
    return _load_path(path, DigestionDocument)


def load_registry_document(path: Path, /) -> RegistryDocument:
    """Load one registry configuration from an explicit JSON path."""
    return _load_path(path, RegistryDocument)


def load_database_build_profile(path: Path, /) -> DatabaseBuildProfileDocument:
    """Load one portable database-build profile from an explicit JSON path."""
    return _load_path(path, DatabaseBuildProfileDocument)


def load_database_build_request(path: Path, /) -> DatabaseBuildRequestDocument:
    """Load one per-run database-build request from an explicit JSON path."""
    return _load_path(path, DatabaseBuildRequestDocument)


def load_protein_input_request(path: Path, /) -> ProteinInputRequestDocument:
    """Load one ordered protein-input preparation request."""
    return _load_path(path, ProteinInputRequestDocument)


def load_derived_protein_input_request(path: Path, /) -> DerivedProteinInputRequestDocument:
    """Load one inventory-to-protein-input derivation request."""
    return _load_path(path, DerivedProteinInputRequestDocument)


def load_decoy_request(path: Path, /) -> DecoyRequestDocument:
    """Load one inventory-to-search-database decoy request."""
    return _load_path(path, DecoyRequestDocument)


def load_decoy_report_request(path: Path, /) -> DecoyReportRequestDocument:
    """Load one peptide-level decoy-method report request."""
    return _load_path(path, DecoyReportRequestDocument)


def load_candidate_request(path: Path, /) -> CandidateRequestDocument:
    """Load one read-only candidate-review request."""
    return _load_path(path, CandidateRequestDocument)


def load_peptide_build_request(path: Path, /) -> PeptideBuildRequestDocument:
    """Load one protein-inventory to peptide-database request."""
    return _load_path(path, PeptideBuildRequestDocument)


def load_peptide_comparison_request(path: Path, /) -> PeptideComparisonRequestDocument:
    """Load one exact peptide-comparison request."""
    return _load_path(path, PeptideComparisonRequestDocument)


def load_uniprot_catalog_request(path: Path, /) -> UniProtCatalogRequestDocument:
    """Load one UniProt catalog request from an explicit JSON path."""
    return _load_path(path, UniProtCatalogRequestDocument)


def load_uniprot_download_request(path: Path, /) -> UniProtDownloadRequestDocument:
    """Load one UniProt download request from an explicit JSON path."""
    return _load_path(path, UniProtDownloadRequestDocument)


def load_builtin_database_build_profile(name: str = "fgcz", /) -> DatabaseBuildProfileDocument:
    """Load one packaged database-build profile by stable name."""
    resource = _PACKAGE_DOCUMENTS.joinpath("build_profiles", name, "profile.json")
    try:
        text = resource.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise ValueError(f"cannot read packaged build profile {name!r}: {error}") from error
    return _load_text(text, str(resource), DatabaseBuildProfileDocument)


def load_builtin_enzyme_document(name: str = "trypsin", /) -> EnzymeDocument:
    """Load one packaged enzyme rule by stable name."""
    resource = _PACKAGE_DOCUMENTS.joinpath("enzymes", name, "rules.json")
    try:
        text = resource.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise ValueError(f"cannot read packaged enzyme document {name!r}: {error}") from error
    return _load_text(text, str(resource), EnzymeDocument)


def load_header_format_document(path: Path, /) -> HeaderFormatDocument:
    """Load one database format from an explicit JSON path."""
    return _load_path(path, HeaderFormatDocument)


def load_header_format_catalog(paths: tuple[Path, ...], /) -> HeaderFormatCatalogDocument:
    """Load an explicit set of independently authored database rules."""
    return HeaderFormatCatalogDocument(
        formats=tuple(load_header_format_document(path) for path in paths)
    )


def load_builtin_header_format_catalog() -> HeaderFormatCatalogDocument:
    """Load every packaged database rule in stable path order."""
    root = _PACKAGE_DOCUMENTS.joinpath("frame_formats")
    rule_resources = sorted(
        (
            child.joinpath("rules.json")
            for child in root.iterdir()
            if child.is_dir() and child.joinpath("rules.json").is_file()
        ),
        key=str,
    )
    return HeaderFormatCatalogDocument(
        formats=tuple(
            _load_text(
                resource.read_text(encoding="utf-8"),
                str(resource),
                HeaderFormatDocument,
            )
            for resource in rule_resources
        )
    )


def _load_path[DocumentT: BaseModel](path: Path, model: type[DocumentT]) -> DocumentT:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise ValueError(f"cannot read JSON document {path}: {error}") from error
    return _load_text(text, str(path), model)


def _load_text[DocumentT: BaseModel](
    text: str,
    source: str,
    model: type[DocumentT],
) -> DocumentT:
    try:
        value = json.loads(text)
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON document {source}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"JSON document {source} must contain an object at its root")
    try:
        return model.model_validate(value)
    except ValidationError as error:
        raise ValueError(f"invalid document {source}: {error}") from error
