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
| Digestion | Enzyme reference, length range, and missed cleavages | Explicit parameter document |
| Source preparation | Ordered target, contaminant, and optional foreign FASTA sources | Explicit parameter document |
| Derived source preparation | Existing source and optional foreign inventories plus role policy | Explicit parameter document |
| Database build profile | Portable naming, metadata, and diagnostic defaults | `documents/build_profiles/fgcz/profile.json` or explicit parameter document |
| Biological build request | Per-run identity, destination, and optional entrapment | Explicit parameter document |
| Effective database build | Fully resolved replay input | Written beside every build |
| Database build result | Artifacts, checksums, summaries, counts, normalization, and generation evidence | Written beside every build |
| Decoy request | Required reverse, shuffle, or DecoyPYrat strategy and output | Explicit parameter document |
| Peptide build/comparison | Digestion, executor, and artifact destinations | Explicit parameter document |
| Candidate review | Comparison threshold, metric, neighbours, and destination | Explicit parameter document |
| UniProt catalog/download | Proteome selection, acquisition mode, destination, and timeout | Explicit parameter document |
| Registry | Backend and indexing/comparison policy | Explicit parameter document |
| Registry diagnostics | Operational FGCZ classifications and decoy prefix | `documents/registry/fgcz.json` |

The packaged diagnostics recognize UniProt, bare UniProt accessions, PDB, RefSeq, GenBank,
Ensembl, and neXtProt identifier shapes. The packaged classifier document contains only broadly
used decoy prefixes (`REV_`, `DECOY_`, `reverse_`) and contaminant prefixes (`CON__`, `CON_`,
`CONTAMINANT_`). Application conventions do not belong here: for example, fasta_gen owns its
`aa|`, `sp|Cont_`, `zh|C...`, and `_p_target` rules in its own JSON document.

Use the packaged documents by omitting options, or consume explicit diagnostic and classifier JSON
beside the command that uses them:

```bash
protein-fasta diagnostics database.fasta \
  --rules rules/diagnostics.json \
  --classifiers rules/classifiers.json
```

The diagnostic document's `schema_version` and `file_version` identify its contract and revision;
`allowed_residues` defines valid normalized sequence symbols; and each
`identifier_namespaces` member gives a reported `name` and matching `pattern`. Each classifier
defines its label/output column plus match and removable prefix/suffix expressions.

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

The CLI consumes the same header-format and classifier documents:

```bash
protein-fasta configured database.fasta proteins.parquet \
  --rules rules/uniprotkb.json \
  --rules rules/refseq.json \
  --classifiers rules/classifiers.json
```

Here each `--rules` document owns one format's detection expression and output columns;
`--classifiers` owns independent labels and removable decorations. The complete validated fields
are described above and in the packaged JSON Schemas.

Committed JSON Schemas for diagnostics, classifiers, header formats, enzymes, source preparation,
biological build, decoys, peptides, candidate review, UniProt, registries, and result evidence are
packaged under `protein_fasta/documents/_schema/`. `make schemas` regenerates them deterministically.
Loading names the invalid source and reports malformed JSON or model-validation errors before
runtime compilation.

## Build configuration

`DatabaseBuildProfileDocument` owns reusable naming, metadata, and diagnostic defaults.
`DatabaseBuildRequestDocument` owns per-run date, naming identity, destination, and optional
biological entrapment.
`resolve_database_build()` applies packaged profile, explicit profile, request, and typed CLI
precedence and returns an `EffectiveDatabaseBuildDocument` with resolved paths. Decoy parameters
are deliberately absent: `DecoyRequestDocument` controls a later operation over the biological
inventory.

The packaged FGCZ profile is ordinary JSON at
`documents/build_profiles/fgcz/profile.json`; copy it when a project needs authored defaults.
Pydantic field defaults remain a final schema fallback, not an invisible application settings
file. The effective request and final `DatabaseBuildResultDocument` are written beside the
biological FASTA. They do not contain decoy policy, registry rows, pair metrics, catalog selection,
installation, or GUI state.

The common build uses the packaged profile and authors its run request from CLI arguments. An
explicit reusable profile is consumed by the same command:

```bash
protein-fasta build protein-input.parquet \
  --output build --project 42261 --dbn 1 --description human
protein-fasta build protein-input.parquet \
  --request build/build.request.json --profile project-profile.json
```

The profile owns naming, metadata, and diagnostic defaults; the generated request owns this run's
output, date, and naming values. The effective document records their resolved combination.

The default metadata grammar constructs the first `aa|<dbname>|...` bookkeeping record with
`CRAPCRAPCRAP` and contaminant section markers with `MRECRAPCRAPCRAP`. These are configuration,
not special cases inside hashing or indexing.

`NamingDocument.allowed_dbname_fields` declares which substitutions its database-name templates
may use; the defaults are `project`, `dbn`, `description`, and `taxid`. Templates may omit any of
them, and absent values collapse through the configured separator. Unknown or attribute/index
expressions, missing filename products, unknown selected templates, and unsupported supplied name
fields are rejected when the Pydantic document loads rather than during a build.

## UniProt configuration

A catalog request selects reference proteomes by default. All proteomes and an explicit provider
query are different storage variants, so a query cannot conflict with an `all` switch. The common
first run authors the request from CLI arguments:

```bash
protein-fasta uniprot-catalog --output catalog
```

It writes this request inside `catalog/` before contacting UniProt:

```json
{
  "schema_version": "0.1",
  "selection": {
    "type": "reference"
  },
  "output_dir": ".",
  "timeout_seconds": 120.0
}
```

`schema_version` selects the request contract; `selection.type` chooses reference proteomes;
`output_dir` is relative to the request file; and `timeout_seconds` records the provider timeout.
Replay it with `protein-fasta uniprot-catalog --request catalog/uniprot-catalog.request.json`.

A download selects exactly one taxon or one proteome identifier and exactly one acquisition mode.
The common proteome-ID form starts directly:

```bash
protein-fasta uniprot-download UP000005640 reviewed --output human-reviewed.fasta
```

Taxonomy selection is an advanced authored request. Save this as `human-reviewed.request.json`:

```json
{
  "schema_version": "0.1",
  "selection": {
    "type": "taxid",
    "taxid": 9606
  },
  "acquisition": {
    "type": "swissprot"
  },
  "output_fasta": "human-reviewed.fasta",
  "timeout_seconds": 120.0
}
```

`selection.type` chooses taxonomy resolution and `taxid` supplies the positive NCBI taxonomy ID;
`acquisition.type` selects reviewed Swiss-Prot entries; `output_fasta` is relative to this request;
and `schema_version` and `timeout_seconds` have the meanings described above. Execute it with:

```bash
protein-fasta uniprot-download --request human-reviewed.request.json
```

Reviewed plus unreviewed changes `acquisition` to `{"type": "swissprot_trembl"}`;
canonical-per-gene uses `{"type": "one_seq_per_gene"}`. The corresponding direct aliases are
`canonical` and `opg`.

Use `{"type": "proteome_id", "proteome_id": "UP000005640"}` instead of the taxid member for
an explicit proteome. Runtime compilation creates one of two resolution behaviors and one of three
acquisition behaviors. Invalid cross-mode field combinations are rejected while loading.

Running `protein-fasta uniprot-download PROTEOME-ID MODE` writes the corresponding authored request
before acquisition; `--request PATH` selects explicit replay without guessing whether a positional
value is a proteome identifier or a file.

## Decoy and peptide configuration

Decoy strategies are separate members. Reverse has no ignored seed or digestion fields; shuffle
adds a seed; DecoyPYrat adds its collision-digestion policy. A peptide request similarly selects
one memory, SQLite, or DuckDB execution member with `workers` and `partition_size`. These passive
documents are compiled once at the workflow root, and every peptide executor returns the same
canonical peptide and mapping schemas.

Digestion policy and enzyme rules remain separate. The first selects an enzyme name, length range,
and missed-cleavage limit; the second gives that named enzyme a versioned `cleavage_pattern`:

```bash
protein-fasta digest database.fasta peptides.parquet \
  --config digestion.json \
  --rules enzymes/trypsin.json
```

In the enzyme JSON, `schema_version` selects the contract, `file_version` identifies the authored
revision, `name` must equal the digestion document's `enzyme`, and `cleavage_pattern` is the regular
expression compiled once before digestion.

The common forms author those documents before running:

```bash
protein-fasta decoy biological-inventory.parquet \
  --output search.fasta --method reverse
protein-fasta peptides search-inventory.parquet \
  --output peptide-products --enzyme trypsin
```

Direct decoy selection covers reverse and shuffle. Use the generated request beside each output as
the starting point for DecoyPYrat, custom seeds, collision-digestion policy, or alternate peptide
executor variants, then replay it with `--request PATH`.

## Registry configuration

`RegistryDocument` contains the creation backend, maximum FASTA size, full-detail entry limit,
metadata-only amino-acid sample size, optional minimum build date, overlap threshold, and the
naming/metadata grammars needed to interpret the collection. The FASTA directory and registry
path are operation arguments rather than hidden deployment defaults.

Registry JSON is reusable policy rather than an automatically authored run request. Supply it
beside the command that consumes it:

```bash
protein-fasta index fasta-directory registry.duckdb --config registry-policy.json
protein-fasta index-inventory search-inventory.parquet registry.duckdb \
  --config registry-policy.json
```

The configured backend controls creation. Reading dispatches from `.sqlite3` or `.duckdb`, so an
existing registry remains readable when a default changes. The registry records its schema,
normalization, hash, fingerprint, and diagnostic-rule versions. Any mismatch that changes stored
meaning requires a full reindex.
