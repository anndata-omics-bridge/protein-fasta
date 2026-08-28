"""Load and compile registry-specific diagnostic rules."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from pydantic import ValidationError

from protein_fasta.compile import make_diagnostic_rules
from protein_fasta.diagnostics.runtime import DiagnosticRules
from protein_fasta.documents import load_builtin_diagnostic_document
from protein_fasta.schema.diagnostics import DiagnosticDocument, EntryClassifierCatalogDocument
from protein_fasta.schema.registry import RegistryDiagnosticDocument

NORMALIZATION_VERSION = "uppercase-one-terminal-stop-v1"
REGISTRY_DIAGNOSTIC_FINGERPRINT_VERSION = "blake2b-128-registry-diagnostics-v1"
_BUILTIN_RESOURCE = resources.files("protein_fasta").joinpath("documents", "registry", "fgcz.json")


@dataclass(frozen=True, slots=True)
class RegistryDiagnosticRules:
    """Schema-free values used by one registry operation."""

    rules: DiagnosticRules
    decoy_prefix: str
    max_reported_id_namespaces: int
    fingerprint: str


def load_registry_diagnostic_document(
    path: Path | None = None,
    /,
) -> RegistryDiagnosticDocument:
    """Load the explicit or packaged registry document."""
    if path is None:
        source = str(_BUILTIN_RESOURCE)
        text = _BUILTIN_RESOURCE.read_text(encoding="utf-8")
    else:
        source = str(path)
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError) as error:
            raise ValueError(f"cannot read registry diagnostics JSON {path}: {error}") from error
    try:
        return RegistryDiagnosticDocument.model_validate_json(text)
    except ValidationError as error:
        raise ValueError(f"invalid registry diagnostics JSON {source}: {error}") from error


def load_registry_diagnostics(path: Path | None = None, /) -> RegistryDiagnosticRules:
    """Compile explicit or packaged registry rules with shared diagnostics."""
    document = load_registry_diagnostic_document(path)

    shared = load_builtin_diagnostic_document()
    combined_diagnostics = DiagnosticDocument(
        file_version=f"{shared.file_version}+registry-{document.file_version}",
        allowed_residues=shared.allowed_residues,
        identifier_namespaces=(*document.identifier_namespaces, *shared.identifier_namespaces),
    )
    classifiers = EntryClassifierCatalogDocument(
        file_version=document.file_version,
        classifiers=document.classifiers,
    )
    payload = json.dumps(
        {
            "normalization_version": NORMALIZATION_VERSION,
            "document": document.model_dump(mode="json"),
            "shared_diagnostics": shared.model_dump(mode="json"),
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    fingerprint = hashlib.blake2b(payload.encode(), digest_size=16).hexdigest()
    return RegistryDiagnosticRules(
        rules=make_diagnostic_rules(combined_diagnostics, classifiers),
        decoy_prefix=document.decoy_prefix,
        max_reported_id_namespaces=document.max_reported_id_namespaces,
        fingerprint=f"blake2b-128:{fingerprint}",
    )
