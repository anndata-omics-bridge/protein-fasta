"""Pure UniProt provider-query construction."""


def reviewed_proteome_query(proteome_id: str, /) -> str:
    """Return the reviewed UniProtKB query for one proteome."""
    return f"(proteome:{proteome_id}) AND (reviewed:true)"


def complete_proteome_query(proteome_id: str, /) -> str:
    """Return the complete UniProtKB query for one proteome."""
    return f"(proteome:{proteome_id})"


def canonical_gene_query(proteome_id: str, /) -> str:
    """Return the GeneCentric query for one proteome."""
    return f"(upid:{proteome_id})"


def reference_taxonomy_query(taxid: int, /) -> str:
    """Return the reference-proteome lookup query for one taxonomy."""
    return f"(taxonomy_id:{taxid}) AND (proteome_type:1)"


def taxonomy_query(taxid: int, /) -> str:
    """Return the fallback proteome lookup query for one taxonomy."""
    return f"(taxonomy_id:{taxid})"
