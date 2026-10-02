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
| `SB` | `SB0140` | SB number: the name of the pattern in the [Mbovis.org](https://www.mbovis.org/) database, or `Not in Mbovis.org` when this database has no name for it (e.g. all *M. tuberculosis* patterns). The spoligotype itself is the `Binary`, `Octal` and `Hexadecimal` codes. With `--db`: the name in that database, or `Not in database` |
| `FileType` | `fastq` | File format: `fastq` (reads) or `fasta` (an assembly, or reads in fasta format: see `Warnings`) |
| `Reads` | `1212727` | Number of reads (or contigs) in the input |
| `Depth` | `65` | Reads only: estimated sequencing depth, all bases divided by 4.4 Mb |
| `MinCount` | `5` | Minimum count used to call a spacer present |
| `Status` | `ok` | `ok`, `warning` (see the next column) or `failed` |
| `Warnings` | | Warnings, separated by `\|`, or the error for failed samples |
| `Species` | `M. bovis` | From the regions of difference and the lineage SNPs, see [Species and lineage](Species-and-lineage) |
| `Lineage` | `BOV` | Most specific lineage of the SNP barcode (e.g. `4.3.4.2`, `2.2.1`, `BOV`), or `mixed: ...` |
| `LineageName` | `M. bovis` | e.g. "Euro-American (LAM)", "East-Asian" (Beijing) |
| `RD9`, `RD4`, `RD1` | `deleted` | `present`, `deleted`, `partial` (some segments of the region missing) or `reduced` (present at a low depth: mixed sample?); also `RD7` and `RD12`, the last columns. See [Species and lineage](Species-and-lineage#species-regions-of-difference) |
| `MTBCFraction` | `1.00` | Reads only: estimated fraction of the reads from the *M. tuberculosis* complex |
| `ClosestSB` | `SB0140 (spacer 7 differs)` | For a pattern not in the database: the closest SB numbers, up to 3 spacers away |
| `RD7`, `RD12` | `deleted` | As `RD9`, `RD4` and `RD1` |
| `SIT` | `SIT451` | Shared international type of the SITVIT2 database: a SIT, `Orphan` (a SITVIT2 pattern without SIT), or `Not in SITVIT2 list` (not among the SITVIT2 patterns of the list, which does not include the SITs created since 2022). Empty without SIT database: see [Installation](Installation#sit-database) |
| `SITVIT2family` | `T-H37Rv` | SITVIT2 spoligotype family of the pattern (also for orphan patterns), e.g. Beijing, LAM3, EAI5, BOV_1 |
| `ClosestSIT` | `SIT451 (spacer 12 differs)` | For a pattern without SIT: the closest SITs, up to 3 spacers away. Not given for a pattern with an SB number (from Mbovis.org, or a `--db` with SB numbers): see [FAQ](FAQ#why-does-my-m-bovis-sample-have-no-sit) |
| `LaLineage` | `La1.8.1` | Lineage of the livestock-associated complex ([Zwyer *et al.* 2021](https://doi.org/10.12688/openreseurope.14029.2)): `La1` (*M. bovis*), `La2` (*M. caprae*), `La3` (*M. orygis*), or an La1 sublineage, e.g. `La1.8.1` (formerly Eu1); `mixed: ...` for a mix; empty for other lineages. See [Species and lineage](Species-and-lineage#livestock-lineages-la1-to-la3) |
| `L1Sublineage` | `L1.2.2.2` | Lineage 1 sublineage ([Netikul *et al.* 2022](https://doi.org/10.1038/s41598-022-05524-0)), in the revised lineage 1 nomenclature, which differs from that of the `Lineage` column (Coll 1.2.1 = L1.2, Coll 1.2.2 = L1.3, Coll 1.1.2 = L1.1.2.2); `mixed: ...` for a mix; empty for other lineages. See [Species and lineage](Species-and-lineage#lineage-1-sublineages) |

Columns are only ever added at the end: the first 6 are those of version 0.2, the next 6 were added in 0.3, the next
8 in 0.4, the next 5 in 0.5, `LaLineage` in 0.8 and `L1Sublineage` in 0.9. The species columns are empty with `--no-species`. Version 0.6 renamed `Spoligotype`
to `SB` and `Closest` to `ClosestSB` (same positions): they hold names in the Mbovis.org database, not the
spoligotype itself.

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
closest patterns, spacers counted from a known variant (`spacer_variants`), species check (with the depth of each region), lineage, livestock lineage and lineage 1 sublineage (`livestock`, `l1`, with
the groups called and the reads supporting each SNP allele) and warnings, and the run information (software versions, parameters, checksums of the reference data).

## MultiQC: `spoligotyping_mqc.json` or `<sample>_spoligotyping_mqc.json`
A [MultiQC custom content](https://docs.seqera.io/multiqc/custom_content) file (JSON, so that octal codes keep their
leading zeros) with the SB number, SIT, octal code,
species, lineage, livestock lineage, lineage 1 sublineage and status of each sample. Run `multiqc` on the output folder to get a "Spoligotyping" section.

## PDF report: `spoligotyping_report.pdf` or `<sample>_spoligotyping.pdf`
Made for quality assurance: everything needed to check a result, and to trace how it was produced.
[Example report](https://github.com/duceppemo/spoligotyper/blob/main/assets/example_report.pdf) (the three samples of the [tutorial](Tutorial)).

![Summary page of the PDF report](https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/report_summary.png)

1. **Summary**: date, operator, and for each sample the octal code, SB number and SIT, species, lineage, pattern and
   status. Samples with the same spoligotype are grouped (the largest groups first), between blue lines, and failed
   samples come last. The warnings and errors follow in a table with a header row, then a
   box for the reviewer's name, date and signature. A table longer than a page continues on the next ones, each part
   titled e.g. "Summary (page 1 of 3)", with its header row repeated.
2. **Samples**: two pages per sample (one for a failed sample, or with `--no-species`), in the order of the summary,
   so their tables are never split (a page too long, with long input paths or many warnings, is scaled down slightly
   to fit):
   * page 1, the spoligotype: octal, hexadecimal and binary codes, and the pattern, its SB number and SIT, the
     closest known patterns when it has no name; the input files: full path (and the real file when it is a symbolic
     link), size, modification date and MD5 checksum; the number of reads and bases, estimated depth, minimum count,
     number of present spacers and their median count; the reads (or contigs) per spacer, present spacers in blue,
     absent spacers seen in some reads in orange; the spacers counted from a known variant (spacer 3 of
     *M. orygis*); and the warnings about the spacers;
   * page 2, the species and lineage: species, lineage with its name and typical spoligotype families, livestock
     lineage, lineage 1 sublineage and amount of MTBC DNA; the regions of difference with their H37Rv coordinates,
     segments found, relative depth and result; the lineage SNPs with their reads; the livestock lineage and lineage
     1 sublineage groups with SNPs carrying the derived allele, with their SNPs covered, with the derived allele and
     with both alleles, their reads and the result (detected, mixed or not detected, explained under each table); and
     the warnings about the species and lineage.
3. **Definitions and methods**: what the codes, SB number, SIT, spoligotype families, lineage, livestock lineage and
   lineage 1 sublineage mean, how the
   regions of difference were determined (H37Rv coordinates, segments and thresholds) and how they compare with PCR
   assays, and how the spacers were counted.
4. **Run information**: operator, user, computer, operating system, start and end time (with time zone), working
   directory, the exact command, parameters, versions of spoligotyper, Python, BBTools and Java, path and MD5
   checksum of the spoligotype database, of the spacer sequences and variants, and of the species and lineage
   reference data, and references.

Every page has the spoligotyper version, the date, user and computer, and "Page x of y" in the footer.

![Sample section of the PDF report](https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/report_sample.png)

MD5 checksums take one to two seconds per GB of input. Use `--no-md5` to skip them (in the PDF and the JSON), and
`--no-pdf` to skip the PDF.
