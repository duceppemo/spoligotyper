# Usage

```
spoligotyper -r1 FILE [-r2 FILE] -o FOLDER [options]
```

## Input files
One sample per run:

| Data | Command |
|---|---|
| Paired-end reads | `spoligotyper -r1 S1_R1.fastq.gz -r2 S1_R2.fastq.gz -o results/` |
| Single-end reads (Illumina, Ion Torrent, nanopore, ...) | `spoligotyper -r1 S1.fastq.gz -o results/` |
| Assembly (contigs or complete genome) | `spoligotyper -r1 S1.fasta -o results/` |

Files can be gzipped or not. The file type is detected from the content, not from the extension.
Use raw or trimmed reads: trimming is not needed.

### Sample name
By default, the sample is named after the `-r1` file, without its extension (`.fastq.gz`, `.fasta`, ...). For fastq
files, the read suffix is removed too: `S1_R1.fastq.gz`, `S1_R1_001.fastq.gz` and `S1_1.fastq.gz` all give `S1`.
Use `-s` to choose another name.

## Options
| Option | Default | Description |
|---|---|---|
| `-r1`, `--r1` | required | Reads (single-end, or R1 of paired-end) or assembly |
| `-r2`, `--r2` | | R2 reads, for paired-end data |
| `-o`, `--output` | required | Folder for the report, created if needed |
| `-s`, `--sample` | from the file name | Sample name, used in the report and its file name |
| `-m`, `--min-count` | 5 for fastq, 1 for fasta | Minimum number of reads containing a spacer to call it present. See [How it works](How-it-works#minimum-count) |
| `-t`, `--threads` | all available | Number of threads |
| `--memory` | `1g` | Java memory for Seal. 1 GB is plenty; see [Troubleshooting](Troubleshooting#java-memory-errors) |
| `--db` | Mbovis.org database | Spoligotype database, see [FAQ](FAQ#can-i-use-another-database) |
| `-v`, `--verbose` | | Also show the Seal command and its output |
| `--version` | | Show the version |

## Output
The report is printed to the screen (standard output) and saved as `<output>/<sample>_spoligotyping.txt`.
Progress and warnings go to standard error, so `spoligotyper ... > result.tsv` saves just the report.
See [Output files](Output-files).

## Many samples
Run spoligotyper once per sample, for example with a loop over paired-end files:
```
for r1 in reads/*_R1.fastq.gz; do
    spoligotyper -r1 "$r1" -r2 "${r1/_R1/_R2}" -o results/ -t 8
done
awk 'FNR == 1 && NR > 1 {next} 1' results/*_spoligotyping.txt > all_samples.tsv
```
Each run takes a few seconds and 1 GB of memory, so several samples can run in parallel (e.g. with GNU parallel).

## Python
spoligotyper can also be used from Python:
```python
from spoligotyper.pipeline import spoligotype

result = spoligotype('S1_R1.fastq.gz', 'S1_R2.fastq.gz', threads=4)
print(result.sample, result.octal, result.spoligotype, result.counts)
```
