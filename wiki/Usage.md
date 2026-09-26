# Usage

```
spoligotyper -r1 FILE [-r2 FILE] -o FOLDER [options]     # One sample
spoligotyper -i FOLDER -o FOLDER [options]               # All the samples in a folder
```

## One sample
| Data | Command |
|---|---|
| Paired-end reads | `spoligotyper -r1 S1_R1.fastq.gz -r2 S1_R2.fastq.gz -o results/` |
| Single-end reads (Illumina, Ion Torrent, ...) | `spoligotyper -r1 S1.fastq.gz -o results/` |
| Assembly (contigs or complete genome). Recommended for nanopore data: see the [FAQ](FAQ#can-i-type-nanopore-reads) | `spoligotyper -r1 S1.fasta -o results/` |

Files can be gzipped or not. The file type is detected from the content, not from the extension. A fasta file with
more than 1,000 sequences and 13 Mb is typed as reads in fasta format rather than as an assembly. A single fastq
file is always typed as single-end reads, even if it contains interleaved pairs.
Use raw or trimmed reads: trimming is not needed.

The sample is named after the `-r1` file, without its extension (`.fastq.gz`, `.fasta`, ...). For fastq files, the
read suffix is removed too: `S1_R1.fastq.gz`, `S1_R1_001.fastq.gz` and `S1_1.fastq.gz` all give `S1`. Use `-s` to
choose another name.

## A folder of samples (batch mode)
```
spoligotyper -i run_2026-09/ -o results/
```
The folder and its subfolders are searched for sequence files (hidden files and folders are skipped, symbolic links
to folders are followed):

| Files | Sample |
|---|---|
| `S1.fasta`, `S1.fa`, `S1.fna`, `S1.fas`, `S1.fsa` (`.gz` or not) | assembly `S1` |
| `S2.fastq`, `S2.fq` (`.gz` or not) | single-end reads `S2` |
| `S3_R1.fastq.gz` + `S3_R2.fastq.gz` (also `_1`/`_2`, `_R1_001`/`_R2_001`) | paired-end reads `S3` |
| `S4_R1.fastq.gz` alone | single-end reads `S4` |

Other files are ignored. If several files give the same sample name (the same file name in two subfolders, `S1.fasta`
and `S1_R1.fastq.gz`, ...), spoligotyper stops and lists them, rather than guessing. The output folder is skipped if
it is inside the input folder.

A sample that fails (e.g. a truncated file) is reported as `failed` in the reports, and the other samples are still
typed; spoligotyper then exits with code 1 so that pipelines notice.

## Options
| Option | Default | Description |
|---|---|---|
| `-r1`, `--r1` | | Reads (single-end, or R1 of paired-end) or assembly |
| `-r2`, `--r2` | | R2 reads, for paired-end data |
| `-i`, `--input` | | Folder of samples (instead of `-r1`) |
| `-s`, `--sample` | from the file name | With `-r1`: sample name |
| `-o`, `--output` | required | Folder for the reports, created if needed |
| `--no-pdf` | | Do not write the PDF report |
| `--no-md5` | | Do not compute the MD5 checksums of the input files (shown in the PDF and JSON reports) |
| `--operator` | user name | Name of the person running the analysis, shown in the PDF report |
| `--no-species` | | Skip the [species check and the lineage](Species-and-lineage): one pass over the reads instead of two |
| `-m`, `--min-count` | 5 for fastq, 1 for fasta | Minimum number of reads containing a spacer to call it present. See [How it works](How-it-works#minimum-count) |
| `--db` | Mbovis.org database | Spoligotype database, see [FAQ](FAQ#can-i-use-another-database) |
| `--sit-db` | downloaded SIT database | SIT database for the `SIT` and `SITVIT2family` columns, see [Installation](Installation#sit-database) |
| `-t`, `--threads` | all available | Number of threads |
| `-j`, `--jobs` | 1 | With `-i`: number of samples typed at the same time, sharing the threads. Each job uses `--memory` |
| `--memory` | `1g` | Java memory for Seal. 1 GB is plenty; see [Troubleshooting](Troubleshooting#java-memory-errors) |
| `-v`, `--verbose` | | Also show the Seal commands and their output |
| `--version` | | Show the version |

## Output
| Mode | Table | PDF report | JSON | MultiQC |
|---|---|---|---|---|
| One sample | `<sample>_spoligotyping.txt` | `<sample>_spoligotyping.pdf` | `<sample>_spoligotyping.json` | `<sample>_spoligotyping_mqc.json` |
| Folder | `spoligotyping.tsv` | `spoligotyping_report.pdf` | `spoligotyping.json` | `spoligotyping_mqc.json` |

The table is also printed to the screen (standard output), and progress and warnings go to standard error, so
`spoligotyper ... > results.tsv` saves just the table. See [Output files](Output-files).

## Python
spoligotyper can also be used from Python:
```python
from spoligotyper.pipeline import spoligotype

result = spoligotype('S1_R1.fastq.gz', 'S1_R2.fastq.gz', threads=4)
print(result.sample, result.octal, result.spoligotype, result.counts, result.warnings)
print(result.species.species, result.lineage.lineage, result.species.mtbc_fraction)
```

## Workflow managers
An [nf-core module](https://github.com/duceppemo/spoligotyper/tree/main/integrations/nf-core) and a
[Galaxy tool](https://github.com/duceppemo/spoligotyper/tree/main/integrations/galaxy) are available in the
repository. The per-sample MultiQC files (`*_mqc.json`) add a "Spoligotyping" table to MultiQC reports: run `multiqc`
on the output folder(s).
