# Architecture

`protein_fasta` owns reusable protein-FASTA acquisition, reading, tabulation, diagnostics,
biological construction, decoy generation, peptide construction, comparison, and registry
persistence. `fasta_gen` is a consumer: it supplies site catalogs and policy, stages or installs
approved artifacts, and presents results in Dash.

The package is an explicit inward dependency structure enforced by Import Linter:

```text
cli
 |
 v
workflow composition roots
(protein_input, database_build, decoy_database, peptide_workflow,
 candidate_analysis, registry_workflow, uniprot_download, uniprot_catalog)
 |
 v
domain and persistence packages
(database, peptide, registry, uniprot, build)
 |
 v
focused computation and I/O
(analytics, reading, validation, diagnostics, frame_formats, schema)
```

This is an allowed-dependency graph, not a required call chain. Root composition modules may know
several inward components; child packages do not import upward or across forbidden siblings. The
exhaustive `.importlinter` contracts are the executable dependency definition.

## Module structure

| Module or package | Responsibility |
| --- | --- |
| `cli.py` | Thin Cyclopts/Loguru delivery adapter |
| `record.py`, `frame.py` | Minimal Python records and optional Polars table products |
| `diagnostic_summary.py`, `summary.py` | Aggregate diagnostic and sequence statistics |
| `reading/`, `validation/` | Lexical I/O, FASTA writing, and fixed normalization |
| `diagnostics/`, `compile.py` | Schema-free diagnostic runtime and document compilation |
| `frame_formats/`, `frame_compile.py` | Polars-native detection and extraction |
| `analytics/` | Backend-free hashes, digestion, comparisons, and clustering |
| `analytics_compile.py` | Enzyme/digestion document compilation |
| `database/`, `inventory.py` | Typed biological/search values and canonical Parquet projection |
| `build/naming.py`, `build/metadata.py` | Database naming and `aa|` metadata construction |
| `build/generation/` | Decoy and biological-entrapment runtime behavior |
| `protein_input.py`, `database_build.py` | Source preparation and biological assembly roots |
| `decoy_database.py`, `decoy_report.py` | Search-database generation and method diagnostics |
| `peptide/`, `peptide_workflow.py` | Peptide model, executors, artifacts, and comparisons |
| `uniprot/`, `uniprot_catalog.py`, `uniprot_download.py` | Provider transport, catalog, resolution, and acquisition |
| `registry/indexing.py` | FASTA/inventory indexing, schema version 11, and typed records |
| `registry/comparisons.py`, `pair_metrics.py` | Storage queries for database comparisons |
| `registry/clustering.py`, `export.py` | Registry-to-analytics projection and stable exports |
| `registry/backend/` | Backend protocol, portable schema, SQLite, DuckDB, and selection |
| `candidate_analysis.py`, `registry_workflow.py` | Read-only review and approved inventory indexing |
| `schema/` | Passive, strict Pydantic storage documents |
| `documents.py`, `documents/` | Explicit-path/resource loading and packaged rules |

## Build stages are explicit

`protein_input.py` prepares ordered sources as canonical Parquet. `database_build.py` consumes
that handoff, renders names and metadata records, normalizes inputs once, rejects conflicting
identifiers, and writes a biological FASTA plus inventory. Optional entrapment belongs to this
biological assembly because its records are members of the biological target space.

`decoy_database.py` subsequently consumes the completed biological inventory and produces the
search FASTA and search inventory. This boundary permits repeatable comparison or replacement of
decoy algorithms without repeating acquisition, source selection, contaminants, or entrapment.

`peptide_workflow.py` consumes either biological or search inventory. Its memory, SQLite, and
DuckDB executors share the same `PeptideExecutor` capability and return the same typed
`PeptideDatabase`; executor-specific persistence never changes the public artifacts.

## Analytics and persistence are separate

`analytics/` creates evidence. It hashes supplied values, digests already-normalized sequences,
calculates similarities, and clusters comparisons. It imports no build code, Pydantic, SQLite,
DuckDB, Polars, Dash, or consuming application. Workflows compile passive documents at their root
and pass exact runtime values inward.

The registry backend protocol expresses the SQL capabilities indexing and comparison queries
exercise. SQLite and DuckDB implement that protocol; engine selection happens once from
configuration during creation and from the file suffix during reading. Canonical inventories can
be indexed directly, while legacy FASTA indexing remains available.

Every full-detail database stores normalized BLAKE2b-128 sequence hashes and materialized target
and contaminant pair statistics. Metadata-only records retain exact aggregate counts without the
entry rows. Registry schema 11 records algorithm/configuration versions and refuses an older
schema; the BLAKE2b cutover requires a full reindex.

## Hash meanings

Hashes are selected by semantic role:

- exact file provenance: MD5, labelled `md5-file-v1` and explicitly non-security;
- protein sequence identity: BLAKE2b-128 over the exact supplied normalized sequence;
- peptide identity: BLAKE2b-128 over the exact supplied peptide;
- ID, description, content, registry-policy, and comparison fingerprints: versioned
  BLAKE2b-128 constructions.

Normalization is never hidden inside a hash function. External source manifests that explicitly
publish SHA-256 remain verified with SHA-256; that authored integrity contract is not redefined.

## Parameters and exchange artifacts

Pydantic models in `schema/` are passive storage documents. They validate serialized parameters
but do not perform FASTA work. Root loading and compilation boundaries construct schema-free
runtime values. The only discriminator branches belong at those construction boundaries;
downstream code invokes the selected behavior.

Each workflow writes effective parameters before computation and typed result evidence after its
data artifacts are durably published. Canonical protein, search, peptide, mapping, comparison,
catalog, and registry exchanges are Parquet; FASTA remains the search-tool exchange; parameter and
result evidence use JSON. Header parsers remain one independently authored document per source
database.

## Optional dependencies and consumers

- Base: records, lexical I/O, diagnostics, summaries, hashing, and schemas.
- `frame`: Polars frame products.
- `uniprot`: HTTP transport for catalog and FASTA acquisition.
- `duckdb`: DuckDB registry and peptide-execution support.
- `generation`: decoy/entrapment generation through `fdr_benchmark`.
- `cli`: Cyclopts, Loguru, Polars/XLSX output, and workflow commands.

The package imports none of `fasta_gen`, Prozor, APB, AnnData, or MuData. Site-specific curated
contaminant/QC catalogs, installation destinations, authorization, and GUI state remain in
`fasta_gen`; reusable UniProt and peptide operations live here.
