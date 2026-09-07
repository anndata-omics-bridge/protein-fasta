# Protein FASTA

`protein-fasta` provides one shared protein-FASTA boundary with deliberately granular APIs:

- a constant-memory Python iterator returning minimal normalized records; and
- an optional Polars frame reader that enriches homogeneous UniProtKB or RefSeq files from
  Pydantic-validated JSON rules;
- backend-free hashing, digestion, comparison, and clustering analytics;
- reproducible source preparation, biological construction, and subsequent decoy generation;
- canonical peptide construction and comparison with memory, SQLite, or DuckDB execution;
- reproducible UniProt proteome catalogs and FASTA acquisition;
- SQLite or optional DuckDB indexing with materialized pair metrics; and
- a short, single-word Cyclopts command surface for each reproducible operation.

**[Online documentation](https://anndata-omics-bridge.github.io/protein-fasta/)**

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

Install `protein-fasta[frame]`. Applications that need one configured database frame select the
accepted public formats when constructing `ProteinDatabase`, then supply the ordered FASTA paths to
`parse()`:

```python
from pathlib import Path

from protein_fasta.frame import ProteinDatabase, refseq, uniprotkb

protein_database = ProteinDatabase(uniprotkb, refseq)
proteins = protein_database.parse(
    (Path("human.fasta"), Path("contaminants.fasta")),
)
```

The result is one row-wise enriched Polars frame. It includes `id`, `description`, `sequence`, the
columns selected by the format profiles, classifier columns, and stable path/checksum/source/record
provenance. Document loading and parser compilation remain private to `protein_fasta`.

Lower-level callers may choose the exact base table or automatic row-wise enrichment:

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
protein-fasta uniprot-download UP000000589 opg
protein-fasta prepare human.fasta protein-input.parquet --id human-uniprot
protein-fasta build protein-input.parquet \
  --output build --project 42261 --dbn 1 --description human
protein-fasta decoy biological.fasta.protein-inventory.parquet \
  --output search.fasta --method reverse
protein-fasta peptides search-inventory.parquet --output peptide-products
protein-fasta prepare --request protein-input.parquet.request.json \
  --output replayed-protein-input.parquet
protein-fasta index databases registry.sqlite3
protein-fasta index-inventory search-inventory.parquet registry.duckdb --config registry.json
protein-fasta registry registry.sqlite3 databases.csv
protein-fasta compare registry.sqlite3 12 comparison.csv
protein-fasta pairs registry.sqlite3 pairs.tsv
protein-fasta cluster registry.sqlite3 clustering.csv
```

Direct workflow commands write their validated authored `*.request.json` before computation.
Replay one explicitly with `--request`; authored requests, resolved `*.effective.json`, observed
`*.result.json`, and data artifacts remain separate.

The CLI also exposes aggregate diagnostics, theoretical digestion, checksums, database builds,
registry indexing, database comparisons, materialized pair exports, and clustering. See the
[build workflows](docs/workflows.md),
[CLI guide and executable walkthrough](docs/cli_walkthrough.md),
[API reference](docs/api.md), and maintained
[architecture](docs/architecture.md).

The package excludes site-specific curated contaminant/QC catalogs, GUI installation workflows,
protein inference, and AnnData/MuData persistence. Those consumers compose the typed artifacts and
APIs without owning FASTA, decoy, peptide, UniProt, or registry computation.

## Development

```bash
uv sync --group dev --extra cli --extra frame --extra duckdb
make check
make docs
```
