# Protein FASTA

`protein-fasta` provides the shared protein-FASTA and reusable database boundary for APB and
`fasta_gen`.

It provides:

- streaming parsing of plain, gzip, and bzip2 FASTA files;
- minimal normalized Python records and compositional diagnostics;
- streaming aggregate diagnostic summaries;
- independent, configurable decoy, contaminant, entrapment, and sentinel labels;
- an exact three-column base Polars frame;
- config-driven homogeneous UniProtKB and RefSeq enrichment;
- versioned fast sequence/peptide hashes and exact-file checksums;
- configured digestion, decoy/entrapment generation, naming, metadata, and database assembly;
- SQLite and DuckDB registry indexing, pair metrics, comparison fingerprints, and clustering; and
- optional table, build, registry, and analytical commands through Cyclopts.

The package deliberately excludes AnnData/MuData persistence, peptide matching, protein
inference, application installation, the curated contaminant/QC catalog, and UniProt download at
this stage.

Continue with the [API reference](api.md), [rule configuration](configuration.md),
[CLI guide](cli.md), [architecture](architecture.md), or [benchmarks](benchmarks.md).
