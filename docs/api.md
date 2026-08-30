# API reference

The API is granular: callers select the smallest record, frame, analytical, build, or registry
operation they need. Creating a `ProteinRecord` never performs database-specific header parsing.

## Installation boundaries

```bash
pip install protein-fasta
pip install 'protein-fasta[frame]'
pip install 'protein-fasta[duckdb]'
pip install 'protein-fasta[generation]'
pip install 'protein-fasta[cli]'
```

## Python and command-line entry points

The package intentionally has no top-level re-export facade. Programmatic callers import the
smallest workflow or computation from its owning module, as in the examples below. The separate
installed command-line entry point is `protein-fasta = protein_fasta.cli:main`; its workflow names
are Cyclopts subcommands rather than additional console scripts.

Request-driven commands accept concise direct arguments for a first run. The CLI validates and
writes the corresponding passive request document before calling the same resolver and runner as a
Python caller. `--request PATH` explicitly replays a saved document, and `--output` may redirect
the replayed products. Python callers construct or load the request document explicitly.

## Records and lexical I/O

```python
from pathlib import Path
from protein_fasta.record import iter_proteins

for protein in iter_proteins(Path("database.fasta.gz")):
    print(protein.id, protein.description, protein.sequence)
```

`ProteinRecord` is exactly `id`, optional `description`, and normalized `sequence`. Formatting
whitespace is removed, letters are upper-cased, and one terminal `*` is removed. The lower-level
`read_records()` and `read_headers()` preserve lexical source meanings when normalization or
sequence materialization is unwanted.

`iter_protein_diagnostics(path, rules)` returns a `ProteinDiagnostics` composition containing the
ordinary record plus raw-header, namespace, overlapping classifications, normalization changes,
and illegal-residue evidence. `summarize_protein_diagnostics(records)` aggregates that stream.

## Polars tables

| Function | Meaning |
| --- | --- |
| `read_basic_protein_frame(path)` | Exact `id`, `description`, `sequence` schema |
| `read_protein_frame(path)` | Row-wise enrichment using built-in formats |
| `read_configured_protein_frame(...)` | Row-wise enrichment using explicit documents |
| `read_strict_protein_frame(path)` | Whole-file built-in enrichment or exact base schema |
| `read_strict_configured_protein_frame(...)` | Whole-file explicit enrichment or exact base schema |
| `read_header_format_diagnostics_frame(...)` | Recognition evidence per candidate format |

Default readers preserve every row and fill parser columns only where exactly one format matches.
Strict readers require one format to match every nonempty row; otherwise they return the exact base
schema.

## Hashing and digestion

```python
from protein_fasta.analytics.hashing import file_checksum, peptide_hash, sequence_hash
from protein_fasta.analytics_compile import make_digestion
from protein_fasta.analytics.digestion import digest_sequence
from protein_fasta.schema.analytics import DigestionDocument

digestion = make_digestion(DigestionDocument(missed_cleavages=1))
peptides = digest_sequence("MPEPTIDEKPEPTIDER", digestion)
```

`sequence_hash()` and `peptide_hash()` return 16-byte BLAKE2b digests of the exact supplied value.
They do not normalize. `file_checksum()` returns hexadecimal MD5 over exact file bytes for
non-security provenance. ID-set, description-set, and ID/sequence content fingerprints are
separately named, versioned BLAKE2b constructions.

Digestion documents select a packaged enzyme rule and peptide-length/missed-cleavage policy.
`digest_sequence()` requires an already-normalized sequence and returns peptide text with its
missed-cleavage count.

## Protein input and biological database construction

Source selection and database assembly are separate APIs. Preparation turns ordered target,
contaminant, and optional foreign FASTA sources into the canonical protein-input Parquet. The
biological build consumes that artifact and never creates decoys.

```python
from pathlib import Path

from protein_fasta.database_build import resolve_database_build, run_database_build
from protein_fasta.documents import (
    load_database_build_profile,
    load_database_build_request,
    load_protein_input_request,
)
from protein_fasta.protein_input import (
    resolve_protein_input_request,
    run_protein_input_preparation,
)

input_request_path = Path("prepare.json")
input_request = resolve_protein_input_request(
    load_protein_input_request(input_request_path),
    request_base=input_request_path.parent,
)
prepared = run_protein_input_preparation(input_request)

profile_path = Path("fgcz.json")
request_path = Path("biological-build.json")
effective = resolve_database_build(
    load_database_build_profile(profile_path),
    load_database_build_request(request_path),
    profile_base=profile_path.parent,
    request_base=request_path.parent,
)
biological = run_database_build(prepared.protein_input_path, effective)
```

Applications such as `fasta_gen` may construct the same Pydantic parameter documents in memory.
Resolvers apply path and precedence rules; runners receive complete effective parameters and
return typed execution values containing the runtime result and durable artifact paths.

`resolve_derived_protein_input_request()` and `run_derived_protein_input_preparation()` prepare an
entrapment source from existing protein or search inventories. They retain target and contaminant
rows, exclude sentinels, prior entrapment, and decoys, and mark retained rows from the optional
foreign inventory with the foreign role.

## Decoy generation

Decoy generation is a subsequent workflow over a completed biological inventory:

```python
from pathlib import Path

from protein_fasta.decoy_database import resolve_decoy_request, run_decoy_generation
from protein_fasta.documents import load_decoy_request

request_path = Path("reverse-decoys.json")
effective = resolve_decoy_request(
    load_decoy_request(request_path),
    request_base=request_path.parent,
)
search = run_decoy_generation(
    Path("build/human.fasta.protein-inventory.parquet"),
    effective,
)
```

The execution returns a `SearchDatabase`, search FASTA, search-inventory Parquet, and typed
effective/result evidence. Reverse, shuffle, and DecoyPYrat are explicit strategy document
variants compiled once by `decoy_compile.make_decoy_generation()` at the workflow boundary.
`DecoyGeneration` is the schema-independent runtime behavior contract shared by database
generation and diagnostic comparison.

The low-level `build_database()` accepts already-resolved target, contaminant, and optional
foreign-source entries plus:

- compiled `DatabaseNaming` for database and filename construction;
- compiled `DatabaseMetadata` for the `aa|` sentinel and contaminant section markers;
- compiled `RegistryDiagnosticRules` for residue checks; and
- optional biological entrapment behavior.

It normalizes once, rejects one identifier with conflicting sequences, deduplicates identical
records, writes deterministic FASTA output, and returns `BiologicalBuildResult`. It does not create
decoys, resolve a contaminant catalog, install into a site collection, or update a registry.

New callers should not assemble that low-level parameter list. Use the workflow API or the
equivalent `build` CLI command. See [Build workflows](workflows.md) for precedence, artifacts,
decoy ownership, and the portrait workflow diagrams.

## Peptide construction and comparison

```python
from pathlib import Path

from protein_fasta.documents import (
    load_peptide_build_request,
    load_peptide_comparison_request,
)
from protein_fasta.peptide_workflow import (
    resolve_peptide_build_request,
    resolve_peptide_comparison_request,
    run_peptide_build,
    run_peptide_comparison,
)

build_request_path = Path("peptides.json")
peptide_request = resolve_peptide_build_request(
    load_peptide_build_request(build_request_path),
    request_base=build_request_path.parent,
)
peptides = run_peptide_build(
    Path("build/human_d.fasta.search-inventory.parquet"),
    peptide_request,
)

comparison_request_path = Path("peptide-comparison.json")
comparison_request = resolve_peptide_comparison_request(
    load_peptide_comparison_request(comparison_request_path),
    request_base=comparison_request_path.parent,
)
comparison = run_peptide_comparison(
    peptides.peptides_path,
    Path("second-peptides.parquet"),
    comparison_request,
)
```

`run_peptide_build()` accepts either protein- or search-inventory Parquet and returns a typed
`PeptideDatabase` plus peptide inventory, protein-peptide mapping, and unique-peptide FASTA.
Memory, SQLite, and DuckDB executors implement the same `PeptideExecutor` capability and produce
the same artifact schemas. Peptide comparison returns one frame with all, target, contaminant,
entrapment, and decoy populations.

## Registry

The registry public surface is split by concern:

| Module | Principal operations |
| --- | --- |
| `registry.indexing` | `rebuild_registry`, `index_fasta`, `index_inventory_entries`, `list_databases`, `get_database` |
| `registry_workflow` | `index_database_inventory` |
| `candidate_analysis` | `resolve_candidate_request`, `run_candidate_analysis` |
| `registry.backend.factory` | `connect`, backend selection by name/suffix |
| `registry.comparisons` | `compare_database`, `compare_candidate`, `find_best_overlap` |
| `registry.pair_metrics` | typed materialized ID/sequence/description/pair counts |
| `registry.export` | `query_similarity_data`, `write_similarity_exports` |
| `registry.clustering` | whole-registry and bounded-neighbourhood clustering |
| `analytics.clustering` | backend-free deterministic average-linkage computation |

`RegistryRecord` exposes counts, normalization evidence, namespace evidence, and versioned target
ID/description/content fingerprints. `DatabaseComparison` exposes directional coverage,
containment, Jaccard metrics, exact pairs, changed shared IDs, and a relationship classification.

Registry schema 11 stores 16-byte sequence hashes and materialized target/contaminant pair rows.
SQLite is the default; DuckDB is selected by configuration when creating and by `.duckdb` suffix
when reading. An old schema is refused and must be fully rebuilt.

`run_candidate_analysis()` compares a canonical inventory against an existing registry without
installing or indexing the candidate. `index_database_inventory()` is the preferred handoff after
approval because it avoids reparsing the generated FASTA.

## UniProt acquisition

`sync_uniprot_catalog()` creates an immutable canonical Parquet catalog from reference, all, or
explicit-query proteome selections. `read_uniprot_catalog()`, `filter_uniprot_catalog()`, and
`latest_uniprot_catalog()` operate locally afterward.

`resolve_uniprot_download()` compiles a taxonomy or proteome-ID selection and a reviewed,
complete, or canonical-gene acquisition request. `run_uniprot_download()` resolves the proteome,
streams and validates the FASTA, publishes it atomically, and returns provider query, release,
reported-count, observed-count, and checksum evidence. Supplying a `UniProtTransport` lets an
application share transport ownership or substitute a test transport.

## Writing and summaries

`write_records(records, path, line_width=60)` serializes lexical `FastaRecord` values.
`summarize_sequences(sequences)` returns exact sequence counts, length statistics, total residues,
and amino-acid frequencies. `SummaryAccumulator` supports streaming and merge operations.
