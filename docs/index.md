# Protein FASTA

`protein-fasta` provides the shared protein-FASTA input boundary for APB and `fasta_gen`.

It provides:

- streaming parsing of plain, gzip, and bzip2 FASTA files;
- generic and UniProt header interpretation;
- configurable target, decoy, contaminant, entrapment, and sentinel classification; and
- protein-sequence normalization, residue validation, and identifier namespaces.

The package deliberately excludes tables, AnnData/MuData persistence, database generation,
digestion, peptide matching, and protein inference.

Continue with the [architecture](architecture.md), or review the initial
[reader benchmarks](benchmarks.md).
