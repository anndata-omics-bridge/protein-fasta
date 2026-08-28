# Rules and frame detection

All maintained regular expressions and configured output columns live in Pydantic-validated JSON
documents. Storage models describe the JSON shape; root composition functions compile them into
schema-free scalar or Polars runtime values.

## Document ownership

| Document | Owns | Packaged location |
| --- | --- | --- |
| Diagnostics | Allowed residues and identifier namespaces | `documents/diagnostics/rules.json` |
| Entry classifiers | Independent labels and removable decorations | `documents/entry_classifiers/rules.json` |
| Header format | Recognition and extracted output columns for one database | `documents/frame_formats/<database>/rules.json` |
| Enzyme | Cleavage expression and stable enzyme identity | `documents/enzymes/<enzyme>/rules.json` |
| Digestion | Enzyme reference, length range, and missed cleavages | Explicit JSON or build document |
| Database build | Inputs, naming, metadata, and generation policy | Explicit JSON |
| Registry | Backend and indexing/comparison policy | Explicit JSON |
| Registry diagnostics | Operational FGCZ classifications and decoy prefix | `documents/registry/fgcz.json` |

The packaged diagnostics recognize UniProt, bare UniProt accessions, PDB, RefSeq, GenBank,
Ensembl, and neXtProt identifier shapes. The packaged classifier document contains only broadly
used decoy prefixes (`REV_`, `DECOY_`, `reverse_`) and contaminant prefixes (`CON__`, `CON_`,
`CONTAMINANT_`). Application conventions do not belong here: for example, fasta_gen owns its
`aa|`, `sp|Cont_`, `zh|C...`, and `_p_target` rules in its own JSON document.

## One rule file per database

Release 0.2 activates exactly two database formats:

| Format | Required enrichment | Optional enrichment |
| --- | --- | --- |
| UniProtKB | `database`, `review_status`, `accession`, `entry_name`, `entry_mnemonic`, `organism_mnemonic`, `protein_name` | `organism_name`, `taxonomy_id`, `gene_name`, `protein_existence`, `sequence_version` |
| RefSeq | `database`, `accession`, `protein_name` | `organism_name` |

The complete UniProt entry name is preserved. `entry_mnemonic` and `organism_mnemonic` are separate
derived columns, while `sp` and `tr` become `reviewed` and `unreviewed` in `review_status`.

Each extracted column declares its name, type, nullability through required/optional placement, and
exactly one source: a literal value or a regular expression with one capture group. A string regex
column may additionally map captured values. Documents cannot replace `id`, `description`, or
`sequence`, and expressions unsupported by Polars are rejected during compilation.

## Detection and decoration handling

Configured classifiers operate on the identifier token. Match expressions add labels without
rewriting it. Removable prefix and suffix expressions both add a label and peel only a temporary
working identifier. The public `id`, `description`, and raw header remain unchanged.

Prefix and suffix rules are applied repeatedly in authored order until a complete pass consumes
nothing. Database detection and extraction then use the fully undecorated working header. This
supports stacked decorations such as `REV_CON__sp|...` while retaining both labels.

The automatic frame reader:

1. reads and normalizes the FASTA once;
2. materializes the internal header column once;
3. applies configured classifiers with native Polars expressions;
4. tests every header against every candidate format; and
5. enriches each row accepted by exactly one format while preserving unmatched rows.

It does not sample or reread the FASTA. Output columns are the stable union of formats that matched
at least one unambiguous row. The strict readers retain the former whole-file contract: one format
must match every row or the result is exactly `id`, `description`, and `sequence`.

## Loading explicit documents

```python
from pathlib import Path

from protein_fasta.documents import (
    load_entry_classifier_document,
    load_header_format_catalog,
)
from protein_fasta.frame import read_configured_protein_frame

catalog = load_header_format_catalog(
    (Path("rules/uniprotkb.json"), Path("rules/refseq.json"))
)
classifiers = load_entry_classifier_document(Path("rules/classifiers.json"))
frame = read_configured_protein_frame(
    Path("database.fasta"),
    catalog,
    classifiers,
)
```

Committed JSON Schemas for diagnostics, classifiers, header formats, enzymes, digestion, database
builds, registries, and registry diagnostics are packaged under
`protein_fasta/documents/_schema/`. `make schemas` regenerates them deterministically. Loading names
the invalid source and reports malformed JSON or model-validation errors before runtime
compilation.

## Build configuration

`DatabaseBuildDocument` separates repeatable file construction from analytics and persistence. It
contains source paths, output directory, build date, naming fields, `NamingDocument`,
`MetadataDocument`, named contaminant blocks, `DecoyDocument`, and optional
`EntrapmentDocument`. It does not contain registry rows, hashes, pair metrics, clustering, catalog
selection, or installation state.

The default metadata grammar constructs the first `aa|<dbname>|...` bookkeeping record with
`CRAPCRAPCRAP` and contaminant section markers with `MRECRAPCRAPCRAP`. These are configuration,
not special cases inside hashing or indexing.

`NamingDocument.allowed_dbname_fields` declares which substitutions its database-name templates
may use; the defaults are `project`, `dbn`, `description`, and `taxid`. Templates may omit any of
them, and absent values collapse through the configured separator. Unknown or attribute/index
expressions, missing filename products, unknown selected templates, and unsupported supplied name
fields are rejected when the Pydantic document loads rather than during a build.

## Registry configuration

`RegistryDocument` contains the creation backend, maximum FASTA size, full-detail entry limit,
metadata-only amino-acid sample size, optional minimum build date, overlap threshold, and the
naming/metadata grammars needed to interpret the collection. The FASTA directory and registry
path are operation arguments rather than hidden deployment defaults.

The configured backend controls creation. Reading dispatches from `.sqlite3` or `.duckdb`, so an
existing registry remains readable when a default changes. The registry records its schema,
normalization, hash, fingerprint, and diagnostic-rule versions. Any mismatch that changes stored
meaning requires a full reindex.
