"""Compile passive biological-database documents into runtime behavior."""

from __future__ import annotations

from typing import TYPE_CHECKING

from protein_fasta.analytics_compile import make_digestion
from protein_fasta.database.metadata import DatabaseMetadata
from protein_fasta.database.naming import DatabaseNaming, DescriptionRules
from protein_fasta.schema.build import (
    BiologicalEntrapmentDocument,
    ForeignSpeciesEntrapmentDocument,
    MetadataDocument,
    NamingDocument,
)

if TYPE_CHECKING:
    from protein_fasta.database.entrapment import EntrapmentGeneration


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
    try:
        from protein_fasta.database.entrapment import (
            make_foreign_species_entrapment_generation,
            make_shuffled_entrapment_generation,
        )
    except ModuleNotFoundError as error:
        missing_name = error.name or ""
        if missing_name == "fdr_benchmark" or missing_name.startswith("fdr_benchmark."):
            raise RuntimeError(
                "entrapment generation requires the 'protein-fasta[generation]' extra"
            ) from error
        raise

    digestion = make_digestion(document.digestion)
    if isinstance(document, ForeignSpeciesEntrapmentDocument):
        return make_foreign_species_entrapment_generation(
            fold=document.fold,
            seed=document.seed,
            enzyme=digestion.cleavage.pattern,
            missed_cleavages=document.digestion.missed_cleavages,
            minimum_length=document.digestion.min_length,
            maximum_length=document.digestion.max_length,
            normalize_i_to_l=document.normalize_i_to_l,
            reject_shared_foreign=document.reject_shared_foreign,
        )
    return make_shuffled_entrapment_generation(
        fold=document.fold,
        seed=document.seed,
        enzyme=digestion.cleavage.pattern,
        missed_cleavages=document.digestion.missed_cleavages,
        minimum_length=document.digestion.min_length,
        maximum_length=document.digestion.max_length,
        normalize_i_to_l=document.normalize_i_to_l,
        fix_peptide_n_term=document.fix_peptide_n_term,
        fix_peptide_c_term=document.fix_peptide_c_term,
    )
