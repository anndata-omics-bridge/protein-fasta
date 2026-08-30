"""Compile passive biological-database documents into runtime behavior."""

from __future__ import annotations

from protein_fasta.analytics_compile import make_digestion
from protein_fasta.database.entrapment import EntrapmentGeneration
from protein_fasta.database.entrapment_generation import (
    make_foreign_species_entrapment_generation,
    make_shuffled_entrapment_generation,
)
from protein_fasta.database.metadata import DatabaseMetadata
from protein_fasta.database.naming import DatabaseNaming, DescriptionRules
from protein_fasta.schema.build import (
    BiologicalEntrapmentDocument,
    ForeignSpeciesEntrapmentDocument,
    MetadataDocument,
    NamingDocument,
)


def make_database_naming(document: NamingDocument, /) -> DatabaseNaming:
    """Compile passive naming configuration into runtime grammar."""
    return DatabaseNaming(
        default_dbname=document.default_dbname,
        dbname=dict(document.dbname),
        filename=dict(document.filename),
        description=DescriptionRules(dict(document.description.mode_flags)),
        decoy_token=document.decoy_token,
        entrapment_token=document.entrapment_token,
        separator=document.separator,
        date_format=document.date_format,
        extension=document.extension,
    )


def make_database_metadata(document: MetadataDocument, /) -> DatabaseMetadata:
    """Compile passive metadata configuration into runtime policy."""
    return DatabaseMetadata(
        id_namespace=document.id_namespace,
        body_sequence=document.body_sequence,
        marker_body_sequence=document.marker_body_sequence,
        date_format=document.date_format,
        provenance_template=document.provenance_template,
        tool=document.tool,
        installer=document.installer,
        org=document.org,
    )


def make_entrapment_generation(
    document: BiologicalEntrapmentDocument,
    /,
) -> EntrapmentGeneration:
    """Compile one entrapment document at the root composition boundary."""
    digestion = make_digestion(document.digestion)
    if isinstance(document, ForeignSpeciesEntrapmentDocument):
        return make_foreign_species_entrapment_generation(
            fold=document.fold,
            seed=document.seed,
            digestion=digestion,
            normalize_i_to_l=document.normalize_i_to_l,
            reject_shared_foreign=document.reject_shared_foreign,
        )
    return make_shuffled_entrapment_generation(
        fold=document.fold,
        seed=document.seed,
        digestion=digestion,
        normalize_i_to_l=document.normalize_i_to_l,
        fix_peptide_n_term=document.fix_peptide_n_term,
        fix_peptide_c_term=document.fix_peptide_c_term,
    )
