"""Regenerate the committed JSON Schemas for authored configuration documents."""

import json
from pathlib import Path

from pydantic import BaseModel

from protein_fasta.schema.analytics import DigestionDocument, EnzymeDocument
from protein_fasta.schema.build import (
    DatabaseBuildProfileDocument,
    DatabaseBuildRequestDocument,
    DatabaseBuildResultDocument,
    EffectiveDatabaseBuildDocument,
)
from protein_fasta.schema.candidate import (
    CandidateRequestDocument,
    CandidateResultDocument,
    EffectiveCandidateRequestDocument,
)
from protein_fasta.schema.decoy import (
    DecoyRequestDocument,
    DecoyResultDocument,
    EffectiveDecoyRequestDocument,
)
from protein_fasta.schema.decoy_report import (
    DecoyReportRequestDocument,
    DecoyReportResultDocument,
    EffectiveDecoyReportDocument,
)
from protein_fasta.schema.diagnostics import DiagnosticDocument, EntryClassifierCatalogDocument
from protein_fasta.schema.frame_formats import HeaderFormatDocument
from protein_fasta.schema.peptide import (
    EffectivePeptideBuildDocument,
    EffectivePeptideComparisonDocument,
    PeptideBuildRequestDocument,
    PeptideBuildResultDocument,
    PeptideComparisonRequestDocument,
    PeptideComparisonResultDocument,
)
from protein_fasta.schema.protein_input import (
    DerivedProteinInputRequestDocument,
    DerivedProteinInputResultDocument,
    ProteinInputRequestDocument,
    ProteinInputResultDocument,
)
from protein_fasta.schema.registry import RegistryDiagnosticDocument, RegistryDocument
from protein_fasta.schema.uniprot import (
    UniProtCatalogRequestDocument,
    UniProtCatalogResultDocument,
    UniProtDownloadRequestDocument,
    UniProtDownloadResultDocument,
)

_SCHEMAS: dict[str, type[BaseModel]] = {
    "candidate_effective.schema.json": EffectiveCandidateRequestDocument,
    "candidate_request.schema.json": CandidateRequestDocument,
    "candidate_result.schema.json": CandidateResultDocument,
    "database_build_effective.schema.json": EffectiveDatabaseBuildDocument,
    "database_build_profile.schema.json": DatabaseBuildProfileDocument,
    "database_build_request.schema.json": DatabaseBuildRequestDocument,
    "database_build_result.schema.json": DatabaseBuildResultDocument,
    "decoy_effective.schema.json": EffectiveDecoyRequestDocument,
    "decoy_request.schema.json": DecoyRequestDocument,
    "decoy_result.schema.json": DecoyResultDocument,
    "derived_protein_input_request.schema.json": DerivedProteinInputRequestDocument,
    "derived_protein_input_result.schema.json": DerivedProteinInputResultDocument,
    "decoy_report_effective.schema.json": EffectiveDecoyReportDocument,
    "decoy_report_request.schema.json": DecoyReportRequestDocument,
    "decoy_report_result.schema.json": DecoyReportResultDocument,
    "diagnostic.schema.json": DiagnosticDocument,
    "digestion.schema.json": DigestionDocument,
    "entry_classifier.schema.json": EntryClassifierCatalogDocument,
    "enzyme.schema.json": EnzymeDocument,
    "header_format.schema.json": HeaderFormatDocument,
    "peptide_build_effective.schema.json": EffectivePeptideBuildDocument,
    "peptide_build_request.schema.json": PeptideBuildRequestDocument,
    "peptide_build_result.schema.json": PeptideBuildResultDocument,
    "peptide_comparison_effective.schema.json": EffectivePeptideComparisonDocument,
    "peptide_comparison_request.schema.json": PeptideComparisonRequestDocument,
    "peptide_comparison_result.schema.json": PeptideComparisonResultDocument,
    "protein_input_request.schema.json": ProteinInputRequestDocument,
    "protein_input_result.schema.json": ProteinInputResultDocument,
    "registry.schema.json": RegistryDocument,
    "registry_diagnostic.schema.json": RegistryDiagnosticDocument,
    "uniprot_catalog_request.schema.json": UniProtCatalogRequestDocument,
    "uniprot_catalog_result.schema.json": UniProtCatalogResultDocument,
    "uniprot_download_request.schema.json": UniProtDownloadRequestDocument,
    "uniprot_download_result.schema.json": UniProtDownloadResultDocument,
}


def main() -> None:
    """Write schemas in deterministic name and JSON-key order."""
    destination = Path("src/protein_fasta/documents/_schema")
    destination.mkdir(parents=True, exist_ok=True)
    for name, model in sorted(_SCHEMAS.items()):
        payload = json.dumps(model.model_json_schema(), indent=2, sort_keys=True) + "\n"
        destination.joinpath(name).write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
