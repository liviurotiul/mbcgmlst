# mbcgmlst

cgMLST allele calling pipeline that maps allele FASTA sequences against assembled genomes and emits a Ridom-style CSV.

## Requirements

- Python 3.9+
- minimap2 in PATH
- For the GUI: a Python build with Tk support (tkinter).

## Basic usage

```bash
mbcgmlst --alleles alleles/ --output results.csv genome.fasta
```

To keep temporary minimap2 files:

```bash
mbcgmlst --alleles alleles/ --output results.csv --keep-temp genome.fasta
```

## GUI

```bash
mbcgmlst
```

## CLI vs GUI

- CLI: run `mbcgmlst` with arguments for automated or batch workflows.
- GUI: run `mbcgmlst-gui`, or just `mbcgmlst` with no arguments to launch the GUI.

### CLI arguments

- `--alleles`: allele FASTA file or directory of FASTA files (one file per locus).
- `--genome-dir`: directory containing genome FASTA files (`.fa`, `.fna`, `.fasta`).
- `--output`: output CSV path (Ridom-style, semicolon-delimited).
- `genomes`: optional list of genome FASTA files (in addition to or instead of `--genome-dir`).
- `--threads`: minimap2 threads (default: 1).
- `--min-identity`: minimum identity for allele hits, 0-1 (default: 0.95).
- `--strict`: require full-length alignments with no indels or clipping.
- `--missing-value`: value for missing loci (default: `NA`).
- `--allow-locus-only`: accept locus-only headers and auto-assign allele numbers.
- `--test`: process only the first genome.
- `--keep-temp`: keep temporary files.

### GUI fields

- Allele FASTA (file or directory)
- Genome FASTA directory
- Output CSV
- Threads (with a Max button)

## Ridom allele compatibility

This pipeline is designed to work with cgMLST allele sets from Ridom (cgmlst.org/ncs).
Common formats are supported:

1) Directory of per-locus FASTA files:
   - File name is the locus (e.g., `ACICU_RS00010.fasta`).
   - Each FASTA header is the allele number only:
     ```
     >26127
     ATGC...
     ```

2) Single combined FASTA:
   - Headers encode both locus and allele:
     ```
     >ACICU_RS00010_26127
     ATGC...
     ```
     or
     ```
     >ACICU_RS00010_allele_26127
     ATGC...
     ```

When using a directory, numeric-only headers are interpreted as allele IDs using the
filename as the locus. The output CSV follows Ridom-style semicolon-delimited format
with locus columns matching the allele set.
