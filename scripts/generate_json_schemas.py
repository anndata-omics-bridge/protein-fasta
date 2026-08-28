"""Regenerate the committed JSON Schemas for authored configuration documents."""

import json
from pathlib import Path

from pydantic import BaseModel

from protein_fasta.schema.analytics import DigestionDocument, EnzymeDocument
from protein_fasta.schema.build import DatabaseBuildDocument
from protein_fasta.schema.diagnostics import DiagnosticDocument, EntryClassifierCatalogDocument
from protein_fasta.schema.frame_formats import HeaderFormatDocument
from protein_fasta.schema.registry import RegistryDiagnosticDocument, RegistryDocument

_SCHEMAS: dict[str, type[BaseModel]] = {
    "database_build.schema.json": DatabaseBuildDocument,
    "diagnostic.schema.json": DiagnosticDocument,
    "digestion.schema.json": DigestionDocument,
    "entry_classifier.schema.json": EntryClassifierCatalogDocument,
    "enzyme.schema.json": EnzymeDocument,
    "header_format.schema.json": HeaderFormatDocument,
    "registry.schema.json": RegistryDocument,
    "registry_diagnostic.schema.json": RegistryDiagnosticDocument,
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
