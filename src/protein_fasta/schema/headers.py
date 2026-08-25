"""Pydantic documents selecting a protein-header convention."""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import Field

from protein_fasta.schema.base import PolicyDocument


class GenericHeaderDocument(PolicyDocument):
    """Use the complete first token as the protein name."""

    kind: Literal["generic"] = "generic"


class UniProtHeaderDocument(PolicyDocument):
    """Interpret UniProt accessions and optional gene names."""

    kind: Literal["uniprot"] = "uniprot"


type HeaderDocument = Annotated[
    GenericHeaderDocument | UniProtHeaderDocument,
    Field(discriminator="kind"),
]
