"""Pydantic documents for configurable registry indexing policy."""

from __future__ import annotations

import datetime
import re
from typing import Literal, Self

from pydantic import Field, model_validator

from protein_fasta.schema.base import DocumentBase
from protein_fasta.schema.build import MetadataDocument, NamingDocument
from protein_fasta.schema.diagnostics import EntryClassifierDocument, NamedPatternDocument


class RegistryBackendDocument(DocumentBase):
    """Storage engine used when a registry is created."""

    backend: Literal["sqlite", "duckdb"] = "sqlite"


class RegistryDocument(DocumentBase):
    """Portable indexing and comparison policy for a FASTA registry.

    Source and destination paths remain operation arguments: they identify one run,
    while this document contains reproducible policy.
    """

    schema_version: Literal["0.1"] = "0.1"
    registry: RegistryBackendDocument = Field(default_factory=RegistryBackendDocument)
    max_fasta_file_size_gib: float = Field(default=5.0, gt=0)
    max_detailed_entries: int = Field(default=200_000, gt=0)
    metadata_aa_sample_size: int = Field(default=20_000, gt=0)
    min_fasta_date: datetime.date | None = None
    overlap_threshold: float = Field(default=0.99, ge=0.0, le=1.0)
    naming: NamingDocument = Field(default_factory=NamingDocument)
    metadata: MetadataDocument = Field(default_factory=MetadataDocument)


class RegistryDiagnosticDocument(DocumentBase):
    """One versioned extension to the shared FASTA diagnostic rules."""

    schema_version: Literal["0.1"] = "0.1"
    file_version: str = Field(min_length=1)
    decoy_prefix: str = Field(min_length=1)
    max_reported_id_namespaces: int = Field(default=32, gt=0)
    identifier_namespaces: tuple[NamedPatternDocument, ...]
    classifiers: tuple[EntryClassifierDocument, ...]

    @model_validator(mode="after")
    def validate_registry_rules(self) -> Self:
        """Require every rule needed by operational registry classification."""
        required_namespaces = frozenset({"fgcz_contaminant", "fgcz_sentinel"})
        namespace_names = {rule.name for rule in self.identifier_namespaces}
        missing_namespaces = required_namespaces - namespace_names
        if missing_namespaces:
            raise ValueError(
                f"missing registry identifier namespaces: {sorted(missing_namespaces)}"
            )

        required_classifiers = frozenset({"contaminant", "decoy", "entrapment", "sentinel"})
        classifier_names = {classifier.name for classifier in self.classifiers}
        missing_classifiers = required_classifiers - classifier_names
        if missing_classifiers:
            raise ValueError(f"missing registry entry classifiers: {sorted(missing_classifiers)}")

        decoy = next(classifier for classifier in self.classifiers if classifier.name == "decoy")
        required_pattern = f"^{re.escape(self.decoy_prefix)}"
        if required_pattern not in decoy.removable_prefix_patterns:
            raise ValueError(
                f"decoy classifier must remove configured decoy_prefix with {required_pattern!r}"
            )
        return self
