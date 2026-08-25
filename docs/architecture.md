# Architecture

`protein_fasta` is the input-format and configured protein-FASTA semantics boundary shared by APB
and `fasta_gen`. Prozor remains independent and accepts ordinary `(protein_id, sequence)` tuples.

```text
consumer application
        |
        +--> protein_fasta.reading          -> FastaRecord
        +--> protein_fasta.compile          -> configured runtime behavior
                    |
                    +--> schema             Pydantic documents
                    +--> headers            plain header values and behavior
                    +--> classification     plain identifier behavior
                    +--> validation         plain sequence behavior
```

The package root is a strict component tree. `compile.py` is the parent composition module; child
packages neither import the parent nor one another. Import Linter enforces that graph.

Pydantic documents represent authored configuration. `compile.py` consumes each discriminator once
and constructs a runtime implementation without a mode field. Runtime methods accept strings,
iterables, compiled regular expressions, and frozen dataclasses. They never receive a Pydantic
model or application container.

FASTA parsing is lexical only. It does not uppercase sequences, strip stops, classify identifiers,
or build application tables. Those are explicit subsequent operations. Database generation,
digestion, registry storage, pandas/Polars tables, AnnData/MuData persistence, peptide matching, and
protein inference belong to consumers.
