# Protein FASTA

`protein-fasta` provides one shared protein-FASTA boundary with deliberately granular APIs:

- a constant-memory Python iterator returning minimal normalized records; and
- an optional Polars frame reader that enriches homogeneous UniProtKB or RefSeq files from
  Pydantic-validated JSON rules;
- backend-free hashing, digestion, comparison, and clustering analytics;
- reproducible database construction with configured naming, metadata, decoys, and entrapment;
- SQLite or optional DuckDB indexing with materialized pair metrics; and
- a short, single-word Cyclopts command surface for each reproducible operation.

The stable high-level Python record is exactly `id`, optional `description`, and normalized
`sequence`:

```python
from pathlib import Path

from protein_fasta.record import iter_proteins

for protein in iter_proteins(Path("proteins.fasta.gz")):
    print(protein.id, protein.sequence)
```

Normalization removes FASTA formatting whitespace, upper-cases the sequence, and removes exactly
one terminal `*`. Callers needing source text use the separately named lower-level
`read_records()` or header-only `read_headers()` operations; raw headers are never reconstructed
from normalized fields.

## Configured diagnostics

`ProteinDiagnostics` composes one `ProteinRecord` with the raw header, identifier namespace,
independent classification labels, normalization changes, and illegal residues:

```python
from pathlib import Path

from protein_fasta.compile import make_diagnostic_rules
from protein_fasta.documents import (
    load_builtin_diagnostic_document,
    load_builtin_entry_classifier_document,
)
from protein_fasta.record import iter_protein_diagnostics

rules = make_diagnostic_rules(
    load_builtin_diagnostic_document(),
    load_builtin_entry_classifier_document(),
)
diagnostics = iter_protein_diagnostics(Path("proteins.fasta"), rules)
```

Classifications overlap: a decoy contaminant can report both labels. Configured decorations are
peeled only from a temporary identifier used for diagnostics and database parsing; the public
record ID remains unchanged. The packaged classifier document contains only format-independent
decorations; application-specific marker conventions belong in an explicit classifier document
supplied to the configured APIs.

For a database-level report, stream those records into the aggregate API:

```python
from protein_fasta.diagnostic_summary import summarize_protein_diagnostics

summary = summarize_protein_diagnostics(diagnostics)
print(summary.namespace_counts, summary.classification_counts)
```

## Polars frames

Install `protein-fasta[frame]`, then choose the exact base table or automatic row-wise enrichment:

```python
from pathlib import Path

from protein_fasta.frame import read_basic_protein_frame, read_protein_frame

base = read_basic_protein_frame(Path("mixed.fasta"))
best = read_protein_frame(Path("uniprot.fasta"))
```

`read_basic_protein_frame()` always returns exactly `id`, `description`, and `sequence`.
`read_protein_frame()` applies the sole matching built-in parser to each row. Unmatched or ambiguous
rows retain their base values and receive null parser fields. `read_strict_protein_frame()` requires
one format to match every row and otherwise returns the exact base schema. Built-in database rules
live in one JSON document per database.

## Table-export CLI

Install the CLI extra and select CSV, TSV, XLSX, or Parquet through the output suffix:

```bash
pip install 'protein-fasta[cli]'
protein-fasta table database.fasta proteins.xlsx --no-sequence
protein-fasta table database.fasta proteins-with-hashes.csv --checksums
protein-fasta table database.fasta proteins.csv --strict
protein-fasta basic database.fasta proteins.csv
protein-fasta formats database.fasta formats.csv
protein-fasta diagnostics database.fasta
protein-fasta digest database.fasta peptides.parquet
protein-fasta checksum database.fasta
protein-fasta build build.json
protein-fasta index databases registry.sqlite3
protein-fasta registry registry.sqlite3 databases.csv
protein-fasta compare registry.sqlite3 12 comparison.csv
protein-fasta pairs registry.sqlite3 pairs.tsv
protein-fasta cluster registry.sqlite3 clustering.csv
```

The CLI also exposes aggregate diagnostics, theoretical digestion, checksums, database builds,
registry indexing, database comparisons, materialized pair exports, and clustering. See the
[build workflows](docs/workflows.md), [CLI guide](docs/cli.md),
[executable CLI walkthrough](docs/cli_walkthrough.md),
[API reference](docs/api.md), and maintained
[architecture](docs/architecture.md).

The package excludes the curated contaminant/QC catalog, GUI installation workflows, UniProt
download at this stage, peptide-to-protein matching, protein inference, and AnnData/MuData
persistence.

## Development

```bash
uv sync --group dev --extra cli --extra frame --extra duckdb --extra generation
make check
make docs
```
