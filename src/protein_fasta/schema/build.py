"""Pydantic documents for database naming and metadata records."""

from __future__ import annotations

import datetime
from enum import StrEnum
from pathlib import Path
from string import Formatter
from typing import Annotated, Literal

from pydantic import Field, model_validator

from protein_fasta.schema.analytics import DigestionDocument
from protein_fasta.schema.artifacts import ArtifactDocument
from protein_fasta.schema.base import DocumentBase


class DescriptionDocument(DocumentBase):
    """Rules for deriving a compact database-description token."""

    mode_flags: dict[str, str] = Field(
        default_factory=lambda: {
            "swissprot": "sp",
            "swissprot_trembl": "",
            "one_seq_per_gene": "1spg",
        }
    )


class NamingDocument(DocumentBase):
    """Grammar for canonical database names and FASTA filenames."""

    default_dbname: str = "project"
    allowed_dbname_fields: tuple[str, ...] = ("project", "dbn", "description", "taxid")
    dbname: dict[str, str] = Field(
        default_factory=lambda: {
            "project": "p{project}_db{dbn}_{description}",
            "fgcz": "fgcz_{taxid}_{description}",
            "derived": "{description}",
        }
    )
    filename: dict[str, str] = Field(
        default_factory=lambda: {
            "decoy": "{dbname}_{decoy_token}_{date}.{extension}",
            "nondecoy": "{dbname}_{date}.{extension}",
            "entrapment_decoy": "{dbname}_{entrapment_token}_{decoy_token}_{date}.{extension}",
            "entrapment_nondecoy": "{dbname}_{entrapment_token}_{date}.{extension}",
        }
    )
    description: DescriptionDocument = Field(default_factory=DescriptionDocument)
    decoy_token: str = "d"
    entrapment_token: str = "e"
    separator: str = "_"
    date_format: str = "%Y%m%d"
    extension: str = "fasta"

    @model_validator(mode="after")
    def validate_templates(self) -> NamingDocument:
        """Reject unresolved template names and unsupported substitutions."""
        if self.default_dbname not in self.dbname:
            raise ValueError(f"default_dbname {self.default_dbname!r} is not a dbname template")
        if len(set(self.allowed_dbname_fields)) != len(self.allowed_dbname_fields):
            raise ValueError("allowed_dbname_fields must not contain duplicates")

        allowed_dbname = set(self.allowed_dbname_fields)
        for name, template in self.dbname.items():
            fields = _template_fields(template, source=f"dbname.{name}")
            unsupported = fields - allowed_dbname
            if unsupported:
                raise ValueError(
                    f"dbname.{name} uses unsupported fields: {sorted(unsupported)}; "
                    f"allowed: {sorted(allowed_dbname)}"
                )

        allowed_filename = {"dbname", "decoy_token", "entrapment_token", "date", "extension"}
        required_filename = {"decoy", "nondecoy", "entrapment_decoy", "entrapment_nondecoy"}
        missing = required_filename - self.filename.keys()
        if missing:
            raise ValueError(f"filename templates are missing: {sorted(missing)}")
        for name, template in self.filename.items():
            fields = _template_fields(template, source=f"filename.{name}")
            unsupported = fields - allowed_filename
            if unsupported:
                raise ValueError(
                    f"filename.{name} uses unsupported fields: {sorted(unsupported)}; "
                    f"allowed: {sorted(allowed_filename)}"
                )
        return self


class MetadataDocument(DocumentBase):
    """Configured database metadata and contaminant-marker record grammar."""

    id_namespace: str = "aa"
    body_sequence: str = "CRAPCRAPCRAP"
    marker_body_sequence: str = "MRECRAPCRAPCRAP"
    date_format: str = "%Y-%m-%d"
    provenance_template: str = "generated w {tool} and installed by {installer}, {org}"
    tool: str = "fasta_gen"
    installer: str = ""
    org: str = "fgcz"


class DecoyMode(StrEnum):
    """Supported protein-decoy generation algorithms."""

    REVERSE = "reverse"
    SHUFFLE = "shuffle"
    DECOYPYRAT = "decoypyrat"


class DecoyDocument(DocumentBase):
    """Validated decoy mode, stochastic seed, and collision digestion."""

    mode: DecoyMode = DecoyMode.REVERSE
    seed: int = 2000
    digestion: DigestionDocument = Field(default_factory=DigestionDocument)


class EntrapmentStrategy(StrEnum):
    """Supported entrapment generation algorithms."""

    SHUFFLED = "shuffled"
    FOREIGN_SPECIES = "foreign-species"


class EntrapmentDocument(DocumentBase):
    """Entrapment strategy, multiplicity, seed, and collision digestion."""

    strategy: EntrapmentStrategy = EntrapmentStrategy.SHUFFLED
    fold: int = Field(default=1, ge=1, le=10)
    seed: int = 2000
    digestion: DigestionDocument = Field(
        default_factory=lambda: DigestionDocument(missed_cleavages=1)
    )
    fix_peptide_n_term: bool = True
    fix_peptide_c_term: bool = True
    normalize_i_to_l: bool = False
    reject_shared_foreign: bool = True


class ShuffledEntrapmentDocument(DocumentBase):
    """Request peptide-shuffled biological entrapment entries."""

    type: Literal["shuffled"] = "shuffled"
    fold: int = Field(default=1, ge=1, le=10)
    seed: int = 2000
    digestion: DigestionDocument = Field(
        default_factory=lambda: DigestionDocument(missed_cleavages=1)
    )
    fix_peptide_n_term: bool = True
    fix_peptide_c_term: bool = True
    normalize_i_to_l: bool = False


class ForeignSpeciesEntrapmentDocument(DocumentBase):
    """Request entrapment proteins selected from prepared foreign-source rows."""

    type: Literal["foreign_species"] = "foreign_species"
    fold: int = Field(default=1, ge=1, le=10)
    seed: int = 2000
    digestion: DigestionDocument = Field(
        default_factory=lambda: DigestionDocument(missed_cleavages=1)
    )
    normalize_i_to_l: bool = False
    reject_shared_foreign: bool = True


BiologicalEntrapmentDocument = Annotated[
    ShuffledEntrapmentDocument | ForeignSpeciesEntrapmentDocument,
    Field(discriminator="type"),
]


class ContaminantBlockDocument(DocumentBase):
    """One named FASTA source inserted as a marked contaminant block."""

    name: str = Field(min_length=1)
    description: str = ""
    path: Path


class DatabaseBuildProfileDocument(DocumentBase):
    """Portable defaults shared by a family of protein database builds."""

    schema_version: Literal["0.3"] = "0.3"
    naming: NamingDocument = Field(default_factory=NamingDocument)
    metadata: MetadataDocument = Field(default_factory=MetadataDocument)
    diagnostics: Path | None = None


class DatabaseBuildRequestDocument(DocumentBase):
    """Per-run sources, identity, destination, and explicit policy overrides."""

    schema_version: Literal["0.3"] = "0.3"
    output_dir: Path
    date: datetime.date
    name_fields: dict[str, str | int]
    template: str | None = None
    naming: NamingDocument | None = None
    metadata: MetadataDocument | None = None
    diagnostics: Path | None = None
    entrapment: BiologicalEntrapmentDocument | None = None
    annotation: str = ""
    installer: str | None = None


class EffectiveDatabaseBuildDocument(DocumentBase):
    """Fully resolved and directly replayable protein database build request."""

    schema_version: Literal["0.3"] = "0.3"
    output_dir: Path
    date: datetime.date
    name_fields: dict[str, str | int]
    template: str
    naming: NamingDocument
    metadata: MetadataDocument
    diagnostics: Path | None = None
    entrapment: BiologicalEntrapmentDocument | None = None
    annotation: str = ""
    installer: str | None = None

    @model_validator(mode="after")
    def validate_naming_request(self) -> EffectiveDatabaseBuildDocument:
        """Reject unresolved naming choices before any source file is read."""
        if self.template not in self.naming.dbname:
            raise ValueError(f"template {self.template!r} is not configured in naming.dbname")
        unsupported = self.name_fields.keys() - set(self.naming.allowed_dbname_fields)
        if unsupported:
            raise ValueError(
                f"name_fields contains unsupported fields: {sorted(unsupported)}; "
                f"allowed: {sorted(self.naming.allowed_dbname_fields)}"
            )
        return self


class DatabaseBuildCountsDocument(DocumentBase):
    """Mutually reconcilable entry counts for one produced FASTA."""

    target: int = Field(ge=0)
    contaminant: int = Field(ge=0)
    entrapment: int = Field(ge=0)
    total: int = Field(ge=0)


class DatabaseBuildNormalizationDocument(DocumentBase):
    """Exact input changes made before database assembly."""

    upper_cased: int = Field(ge=0)
    terminal_stops_stripped: int = Field(ge=0)
    duplicates_dropped: int = Field(ge=0)


class DatabaseBuildSummaryDocument(DocumentBase):
    """Sequence-length and amino-acid summary persisted with a build."""

    n_sequences: int = Field(ge=0)
    length_min: int | None = Field(default=None, ge=0)
    length_max: int | None = Field(default=None, ge=0)
    length_mean: float | None = Field(default=None, ge=0)
    length_q1: float | None = Field(default=None, ge=0)
    length_median: float | None = Field(default=None, ge=0)
    length_q3: float | None = Field(default=None, ge=0)
    total_residues: int = Field(ge=0)
    aa_counts: dict[str, int]
    aa_frequencies: dict[str, float]


class DatabaseBuildEntrapmentEvidenceDocument(DocumentBase):
    """Entrapment strategy identity and achieved multiplicity."""

    strategy: str
    seed: int
    requested_fold: int = Field(ge=1)
    achieved_fold: int = Field(ge=0)
    failures: int = Field(ge=0)
    proteins_affected: int = Field(ge=0)
    source_proteins: int = Field(ge=0)


class DatabaseBuildResultDocument(DocumentBase):
    """Versioned machine-readable evidence for one completed database build."""

    schema_version: Literal["0.3"] = "0.3"
    protein_fasta_version: str
    effective_request: EffectiveDatabaseBuildDocument
    input_artifact: ArtifactDocument
    biological_fasta: ArtifactDocument
    protein_inventory: ArtifactDocument
    sidecar_artifacts: tuple[ArtifactDocument, ...] = ()
    counts: DatabaseBuildCountsDocument
    normalization: DatabaseBuildNormalizationDocument
    summary: DatabaseBuildSummaryDocument
    entrapment: DatabaseBuildEntrapmentEvidenceDocument | None = None


def _template_fields(template: str, *, source: str) -> set[str]:
    """Return plain replacement fields from one authored format template."""
    try:
        parsed = Formatter().parse(template)
        fields = {field for _, field, _, _ in parsed if field is not None}
    except ValueError as error:
        raise ValueError(f"{source} is not a valid format template: {error}") from error
    unsafe = {field for field in fields if not field.isidentifier()}
    if unsafe:
        raise ValueError(f"{source} uses unsafe fields: {sorted(unsafe)}")
    return fields
