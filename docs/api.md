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

## Database construction

The stable workflow API is `resolve_database_build()` followed by `run_database_build()`. The
first call combines a `DatabaseBuildProfileDocument`, a `DatabaseBuildRequestDocument`, their
base directories, and optional typed `DatabaseBuildOverrides`. The second consumes the complete
`EffectiveDatabaseBuildDocument` and returns `DatabaseBuildExecution`.

```python
from pathlib import Path

from protein_fasta.database_build import resolve_database_build, run_database_build
from protein_fasta.documents import (
    load_database_build_profile,
    load_database_build_request,
)

profile_path = Path("fgcz.json")
request_path = Path("request.json")
effective = resolve_database_build(
    load_database_build_profile(profile_path),
    load_database_build_request(request_path),
    profile_base=profile_path.parent,
    request_base=request_path.parent,
)
execution = run_database_build(effective)
```

This is the API `fasta_gen` should call. It may construct the Pydantic documents in memory instead
of loading JSON; the resolver and execution call stay identical. The returned execution contains
the runtime `PipelineResult`, typed `DatabaseBuildResultDocument`, and both JSON paths.

The internal low-level `build_database()` currently accepts already-resolved target,
contaminant, and optional foreign-source entries plus:

- `NamingDocument` for database and filename construction;
- `MetadataDocument` for the `aa|` sentinel and contaminant section markers;
- compiled `RegistryDiagnosticRules` for residue checks and the decoy prefix;
- `DecoyDocument` and optional `EntrapmentDocument` for generation behavior.

It normalizes once, rejects one identifier with conflicting sequences, deduplicates identical
records, generates requested records, writes deterministic FASTA output, and returns
`PipelineResult`. It does not resolve a contaminant catalog, install into a site collection, or
update a registry; those are application workflows.

New callers should not assemble that low-level parameter list. Use the workflow API or the
equivalent `build` CLI command. See [Build workflows](workflows.md) for precedence, artifacts,
decoy ownership, and the portrait workflow diagrams.

## Registry

The registry public surface is split by concern:

| Module | Principal operations |
| --- | --- |
| `registry.indexing` | `rebuild_registry`, `index_fasta`, `list_databases`, `get_database` |
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

## Writing and summaries

`write_records(records, path, line_width=60)` serializes lexical `FastaRecord` values.
`summarize_sequences(sequences)` returns exact sequence counts, length statistics, total residues,
and amino-acid frequencies. `SummaryAccumulator` supports streaming and merge operations.
