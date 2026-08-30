# Architecture

`protein_fasta` owns reusable protein-FASTA acquisition, reading, tabulation, diagnostics,
biological construction, decoy generation, peptide construction, comparison, and registry
persistence. `fasta_gen` is a consumer: it supplies site catalogs and policy, stages or installs
approved artifacts, and presents results in Dash.

## Delivery and public entry points

The installed distribution declares one console script:

```toml
[project.scripts]
protein-fasta = "protein_fasta.cli:main"
```

Cyclopts exposes the operations beneath that entry point as subcommands; `build`, `decoy`,
`prepare`, `peptides`, `uniprot-download`, and the registry commands are not separate executables.
Programmatic callers import the workflow functions from their owning modules. The package has no
top-level re-export facade: empty package initializers keep ownership and dependency direction
visible.

The package has an explicit inward dependency structure enforced by Import Linter:

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
(database, peptide, registry, uniprot)
 |
 v
focused computation and I/O
(analytics, reading, validation, diagnostics, frame_formats, schema)
```

This is an allowed-dependency graph, not a required call chain. Root composition modules may know
several inward components. The exhaustive `.importlinter` contracts enforce the current layered
DAG and the framework/consumer isolation boundaries; they do not yet prove the stronger directed-
folder rule described below.

## Module structure

| Module or package | Responsibility |
| --- | --- |
| `cli.py` | Single Cyclopts/Loguru delivery and composition adapter |
| `artifact_io.py` | Atomic artifact publication, JSON persistence, and artifact evidence |
| `record.py`, `frame.py` | Minimal Python records and optional Polars table products |
| `diagnostic_summary.py`, `summary.py` | Aggregate diagnostic and sequence statistics |
| `reading/`, `validation/` | Lexical I/O, FASTA writing, and fixed normalization |
| `diagnostics/`, `compile.py` | Schema-free diagnostic runtime and document compilation |
| `frame_formats/`, `frame_compile.py` | Polars-native detection and extraction |
| `analytics/` | Backend-free hashes, digestion, comparisons, and clustering |
| `analytics_compile.py` | Enzyme/digestion document compilation |
| `database/`, `inventory.py` | Typed biological/search values and canonical Parquet projection |
| `database/naming.py`, `database/metadata.py` | Database naming and `aa|` metadata construction |
| `database/decoy*.py`, `database/entrapment.py` | Decoy and biological-entrapment runtime behavior |
| `database_compile.py`, `decoy_compile.py` | Passive build/strategy document compilation |
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

## Current dependency-structure status

The current graph is acyclic and all Import Linter contracts pass. It is not yet a strict directed
folder tree in which every child imports no parent module and at most one sibling package. The
remaining structural edges are explicit:

| Importing child | Current dependency | Directed-folder consequence |
| --- | --- | --- |
| `registry/` | parent-root `compile`, `documents`, `record`, and `summary` | upward imports remain |
| `peptide/` | sibling `analytics/` and sibling `registry/` | more than one direct sibling |
| `registry/` | sibling `analytics/`, `database/`, `diagnostics/`, `reading/`, and `schema/` | more than one direct sibling |

`uniprot/ -> reading/` uses one directed sibling edge and therefore does not violate the
sibling-count rule. The existing architecture tests additionally keep every
initializer empty, name the root modules allowed to compose multiple children, and constrain the
CLI's imported package components. A future folder migration must move ownership or composition;
adding forwarding modules solely to hide these edges would not improve the architecture.

## Runtime behavior and storage documents

The UniProt and peptide workflows compile passive discriminator documents once at their root into
runtime objects that own resolution, acquisition, or execution behavior. Those branches are
storage-boundary factories, not missing methods on the Pydantic documents.

`decoy_compile.py` is the single storage-to-runtime composition owner for reverse, shuffle, and
DecoyPYrat documents. It constructs the schema-independent `DecoyGeneration` behavior used by
both database generation and diagnostic comparison. `database_compile.py` performs the equivalent
one-time compilation for naming, metadata, and biological entrapment documents. The child
`database/` package therefore has no dependency on Pydantic schemas or workflow roots.

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
but do not perform FASTA work. Root loading and compilation boundaries construct runtime values.
Discriminator branches belong at those construction boundaries; downstream code invokes the
selected behavior. The duplicate decoy compilation bridge documented above is the remaining
exception, not the intended pattern.

Each workflow writes effective parameters before computation and typed result evidence after its
data artifacts are durably published. Canonical protein, search, peptide, mapping, comparison,
catalog, and registry exchanges are Parquet; FASTA remains the search-tool exchange; parameter and
result evidence use JSON. Header parsers remain one independently authored document per source
database.

The CLI normally starts from concise direct arguments, constructs the same passive request document,
and persists it without replacement before calling the shared resolver and runner. `--request`
explicitly selects replay; an output override may redirect products, while direct scientific options
cannot be mixed into replay. This keeps interactive use reproducible without introducing a second
workflow API or inferring whether a positional string is a request path.

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
