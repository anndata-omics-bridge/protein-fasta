# Architecture

`protein_fasta` owns reusable protein-FASTA reading, tabulation, diagnostics, analytics,
database construction, and registry persistence. `fasta_gen` is a consumer: it selects curated
catalog inputs, stages and installs builds, and presents results in Dash.

The package is an explicit inward dependency structure enforced by Import Linter:

```text
cli
 |
 v
database_build | diagnostic_summary
 |
 v
registry
 |
 v
build | record | frame | writing | summary
 |
 v
analytics_compile
 |
 v
compile | frame_compile | documents
 |
 v
analytics | reading | validation | diagnostics | frame_formats | schema
```

This is an allowed-dependency graph, not a required call chain. Modules on the same line are
independent. Root composition modules may know several inward components; child packages do not
import upward or sideways.

## Module structure

| Module or package | Responsibility |
| --- | --- |
| `cli.py` | Thin Cyclopts/Loguru delivery adapter |
| `record.py`, `frame.py` | Minimal Python records and optional Polars table products |
| `diagnostic_summary.py`, `summary.py` | Aggregate diagnostic and sequence statistics |
| `reading/`, `validation/`, `writing.py` | Lexical I/O and fixed normalization |
| `diagnostics/`, `compile.py` | Schema-free diagnostic runtime and document compilation |
| `frame_formats/`, `frame_compile.py` | Polars-native detection and extraction |
| `analytics/` | Backend-free hashes, digestion, comparisons, and clustering |
| `analytics_compile.py` | Enzyme/digestion document compilation |
| `build/naming.py` | Database-name and filename construction only |
| `build/metadata.py` | `aa|` metadata and contaminant-marker construction only |
| `build/generation/` | Decoy and entrapment runtime behavior |
| `database_build.py` | Public assembly composition root |
| `registry/indexing.py` | FASTA scanning, indexing, schema version 11, and typed records |
| `registry/comparisons.py`, `pair_metrics.py` | Storage queries for database comparisons |
| `registry/clustering.py`, `export.py` | Registry-to-analytics projection and stable exports |
| `registry/backend/` | Backend protocol, portable schema, SQLite, DuckDB, and selection |
| `registry/filenames.py`, `metadata.py`, `classification.py` | Existing-file interpretation and operational kinds |
| `schema/` | Passive, strict Pydantic JSON documents |
| `documents.py`, `documents/` | Explicit-path/resource loading and packaged rules |

## Build and analytics are separate

`build/` creates database artifacts. It renders names, metadata records, contaminant block
markers, decoys, and entrapment records. `database_build.py` normalizes inputs once, rejects
conflicting identifiers, assembles records deterministically, writes the FASTA, and returns a
typed result.

`analytics/` creates evidence. It hashes supplied values, digests already-normalized sequences,
calculates similarities, and clusters comparisons. It imports no build code, Pydantic, SQLite,
DuckDB, Polars, Dash, or consuming application. A build may consume an inward digestion runtime;
analytics never imports build behavior.

This separation is why filename construction and `aa|... CRAPCRAP...` bookkeeping are not mixed
with sequence hashes or comparison fingerprints.

## Registry boundary

The registry is reusable infrastructure, not GUI state. Its backend protocol expresses the SQL
capabilities indexing and comparison queries exercise. SQLite and DuckDB implement that protocol;
engine selection happens once from configuration during creation and from the file suffix during
reading.

Every full-detail database stores normalized BLAKE2b-128 sequence hashes and materialized target
and contaminant pair statistics. Metadata-only records retain exact aggregate counts without the
entry rows. Comparison calculations and average-linkage clustering consume typed values outside
the backend. Registry schema 11 records algorithm/configuration versions and refuses an older
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

## Configuration and runtime behavior

Pydantic models in `schema/` are storage documents. They validate JSON but do not perform FASTA
work. Root loading/compilation boundaries construct schema-free runtime values. The only
discriminator branches belong at those construction boundaries; downstream code invokes the
selected behavior.

Database-build profile and request JSON are resolved once into a complete effective document.
`run_database_build()` writes that document before computation and a typed result afterward.
Registry JSON contains backend and indexing policy while source and destination paths remain
operation arguments. Header parsers remain one independently authored JSON file per source database.

## Optional dependencies and consumers

- Base: records, lexical I/O, diagnostics, summaries, hashing, and schemas.
- `frame`: Polars frame products.
- `duckdb`: DuckDB registry adapter and Arrow transfer support.
- `generation`: decoy/entrapment generation through `fdr_benchmark`.
- `cli`: Cyclopts, Loguru, Polars/XLSX output, and generation commands.

The package imports none of `fasta_gen`, Prozor, APB, AnnData, or MuData. UniProt download and the
curated contaminant/QC catalog remain in `fasta_gen` at this stage; neither is part of the registry
or analytical core described here.
