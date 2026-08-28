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
- UniProt catalog synchronization and reproducible proteome FASTA acquisition;
- separate source preparation, biological assembly, and subsequent decoy generation;
- canonical peptide inventories, mappings, FASTA products, and exact comparisons;
- SQLite and DuckDB registry indexing, pair metrics, comparison fingerprints, and clustering; and
- optional table, build, registry, and analytical commands through Cyclopts.

The package deliberately excludes AnnData/MuData persistence, protein inference, application
installation and authorization, and the site-curated contaminant/QC catalog.

Start with the [build workflows](workflows.md), then continue with the [API reference](api.md),
[artifact contracts](artifacts.md), [rule configuration](configuration.md), [CLI guide](cli.md),
[architecture](architecture.md), or [benchmarks](benchmarks.md).
