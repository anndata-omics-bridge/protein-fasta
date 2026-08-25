# Protein FASTA

`protein-fasta` supplies the protein-FASTA input boundary shared by APB and `fasta_gen`:

- streaming plain, gzip, and bzip2 record parsing;
- generic and UniProt header interpretation;
- configured target/decoy/contaminant/entrapment/sentinel classification; and
- protein-sequence normalization, residue validation, and identifier namespaces.

The import package is `protein_fasta`. It deliberately contains no pandas, Polars, NumPy, AnnData,
MuData, database-building, digestion, or protein-inference code. See
[`docs/architecture.md`](docs/architecture.md).

## Reading records

```python
from pathlib import Path

from protein_fasta.headers.generic import parse_header
from protein_fasta.reading.parser import read_records

for record in read_records(Path("proteins.fasta.gz")):
    header = parse_header(record.raw_header)
    print(header.identifier, record.sequence)
```

Use `parse_records(lines)` for an already-open text stream and `parse_text(text)` for explicit
inline FASTA content. These separate functions keep ownership of the stream and the meaning of a
string unambiguous.

## Compiling configured behavior

```python
from protein_fasta.compile import make_header_interpreter, make_sequence_validator
from protein_fasta.schema.headers import UniProtHeaderDocument
from protein_fasta.schema.validation import SequenceValidationDocument

header_interpreter = make_header_interpreter(UniProtHeaderDocument())
sequence_validator = make_sequence_validator(SequenceValidationDocument())

header = header_interpreter.interpret("sp|P12345|KINASE Protein kinase GN=MAPK1")
sequence = sequence_validator.normalize("mpeptidek*")
```

The Pydantic documents stop at the `make_*` boundary. The resulting runtime objects retain only
the values needed by their operations.

See the maintained [architecture](docs/architecture.md) and initial [reader benchmark](docs/benchmarks.md).

## Development

```bash
uv sync --group dev
make check
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
```

All Python commands run from the synchronized project `.venv`.
