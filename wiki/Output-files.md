# Output files

## Table: `spoligotyping.tsv` or `<sample>_spoligotyping.txt`
A tab-separated table with a header and one line per sample. The same table is printed to the screen.

| Column | Example | Description |
|---|---|---|
| `Sample` | `AF2122_97` | Sample name |
| `SpacerCount` | `56:47:0:58:...` | Number of reads containing each spacer, spacers 1 to 43, separated by `:`. For assemblies, the number of contigs |
| `Binary` | `1101101000001...` | 43 digits: 1 = spacer present (count ≥ `MinCount`), 0 = absent |
| `Octal` | `664073777777600` | 15-digit octal code |
| `Hexadecimal` | `6D-03-5F-7F-FF-60` | Hexadecimal code, 6 blocks |
| `Spoligotype` | `SB0140` | SB number of the pattern in the [Mbovis.org](https://www.mbovis.org/) database, or `Spoligo not found` |
| `FileType` | `fastq` | File format: `fastq` (reads) or `fasta` (an assembly, or reads in fasta format: see `Warnings`) |
| `Reads` | `1212727` | Number of reads (or contigs) in the input |
| `Depth` | `65` | Reads only: estimated sequencing depth, all bases divided by 4.4 Mb |
| `MinCount` | `5` | Minimum count used to call a spacer present |
| `Status` | `ok` | `ok`, `warning` (see the next column) or `failed` |
| `Warnings` | | Warnings, separated by `\|`, or the error for failed samples |
| `Species` | `M. bovis` | From the regions of difference and the lineage SNPs, see [Species and lineage](Species-and-lineage) |
| `Lineage` | `BOV` | Most specific lineage of the SNP barcode (e.g. `4.3.4.2`, `2.2.1`, `BOV`), or `mixed: ...` |
| `LineageName` | `M. bovis` | e.g. "Euro-American (LAM)", "East-Asian" (Beijing) |
| `RD9`, `RD4`, `RD1` | `deleted` | `present`, `deleted` or `partial` (also `RD7` and `RD12`, the last two columns) |
| `MTBCFraction` | `1.00` | Reads only: estimated fraction of the reads from the *M. tuberculosis* complex |
| `Closest` | `SB0140 (spacer 7 differs)` | For a pattern not in the database: the closest SB numbers, up to 3 spacers away |
| `RD7`, `RD12` | `deleted` | `present`, `deleted` or `partial` |
| `SIT` | `SIT451` | Shared international type of the SITVIT2 database: a SIT, `Orphan` (a SITVIT2 pattern without SIT), or `Spoligo not found`. Empty without SIT database: see [Installation](Installation#sit-database) |
| `SITVIT2family` | `T-H37Rv` | SITVIT2 spoligotype family of the pattern (also for orphan patterns), e.g. Beijing, LAM3, EAI5, BOV_1 |
| `ClosestSIT` | `SIT451 (spacer 12 differs)` | For a pattern without SIT: the closest SITs, up to 3 spacers away |

Columns are only ever added at the end: the first 6 are those of version 0.2, the next 6 were added in 0.3, the next
8 in 0.4 and the last 5 in 0.5. The species columns are empty with `--no-species`.

### Octal code
The binary pattern is cut into 14 groups of 3 spacers, and each group is written as one octal digit (000 = 0,
001 = 1, ..., 111 = 7). Spacer 43 is the 15th digit, 0 or 1. It is the standard code used by SITVIT and in
publications, and it can be converted back to the binary pattern.

### Hexadecimal code
The binary pattern is cut into 6 blocks of 7, 7, 7, 7, 8 and 7 spacers, and each block is written as a 2-digit
hexadecimal number.

### Reading the counts
`SpacerCount` shows how confident each call is:
* **Reads**: present spacers usually have tens of reads and absent spacers 0. Counts close to `MinCount`
  deserve a second look, and spoligotyper warns about absent spacers seen in a few reads. See
  [How it works](How-it-works#minimum-count).
* **Paired-end reads**: when a spacer is found in one read of a pair, both reads are counted, so counts are about
  twice the number of DNA fragments containing the spacer.
* **Assemblies**: counts are 0 or 1 (occasionally 2, if a spacer is split over two contigs or repeated).

## JSON: `spoligotyping.json` or `<sample>_spoligotyping.json`
Everything in the table and the PDF report, for pipelines: for each sample the input files, codes, spacer counts,
closest patterns, species check (with the depth of each region), lineage (with the reads supporting each SNP
allele) and warnings, and the run information (software versions, parameters, checksums of the reference data).

## MultiQC: `spoligotyping_mqc.json` or `<sample>_spoligotyping_mqc.json`
A [MultiQC custom content](https://docs.seqera.io/multiqc/custom_content) file (JSON, so that octal codes keep their
leading zeros) with the spoligotype, octal code,
species, lineage and status of each sample. Run `multiqc` on the output folder to get a "Spoligotyping" section.

## PDF report: `spoligotyping_report.pdf` or `<sample>_spoligotyping.pdf`
Made for quality assurance: everything needed to check a result, and to trace how it was produced.
[Example report](https://github.com/duceppemo/spoligotyper/blob/main/assets/example_report.pdf) (the three samples of the [tutorial](Tutorial)).

![Summary page of the PDF report](https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/report_summary.png)

1. **Summary**: date, operator, and for each sample the spoligotype, octal code, species, lineage, pattern and
   status. Samples with the same spoligotype are grouped (the largest groups first; alternate groups shaded), and
   failed samples come last. Warnings and errors are listed below, followed by a box for the reviewer's name, date
   and signature.
2. **Samples**: one page per sample, in the order of the summary, so its tables are never split, with
   * the spoligotype, octal, hexadecimal and binary codes, and the pattern;
   * the input files: full path (and the real file when it is a symbolic link), size, modification date and MD5
     checksum;
   * the number of reads and bases, estimated depth, minimum count, number of present spacers and their median count;
   * the closest known patterns, when the pattern is not in the database;
   * the species and lineage: regions of difference with their relative depth, lineage with its name and typical
     spoligotype families, amount of MTBC DNA, and the reads supporting each lineage SNP;
   * the reads (or contigs) per spacer: present spacers in blue, absent spacers seen in some reads in orange;
   * the warnings, or the error of a failed sample.
3. **Run information**: operator, user, computer, operating system, start and end time (with time zone), working
   directory, the exact command, parameters, versions of spoligotyper, Python, BBTools and Java, path and MD5
   checksum of the spoligotype database and of the spacer sequences, the method, and references.

Every page has the spoligotyper version, the date, user and computer, and "Page x of y" in the footer.

![Sample section of the PDF report](https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/report_sample.png)

MD5 checksums take one to two seconds per GB of input. Use `--no-md5` to skip them (in the PDF and the JSON), and
`--no-pdf` to skip the PDF.
