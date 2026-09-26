# Changelog

## 0.5.0 (2026-09-25)
* Species check with five regions of difference, as in the classical RD PCR scheme: RD7 and RD12 are added to RD1,
  RD4 and RD9. New species calls: *M. africanum* lineage 5 or 6, *M. orygis* or *M. caprae*, *M. microti*,
  *M. pinnipedii* or *M. mungi*, *M. canettii* (also from the absence of standard spacers), and the Dassie bacillus.
  See [Species and lineage](Species-and-lineage).
* Fixed: an RD1 deletion with RD4 present was reported as "e.g. *M. microti*": *M. microti* is RD1 present, the
  Dassie bacillus is RD1 deleted.
* Fixed: *M. canettii* strains with the lineage 4 allele were called *M. tuberculosis*, and those without RD4 could
  have been called *M. bovis*.
* No "RD4 present" warning for *M. caprae* and *M. orygis*, which carry the *M. bovis* clade SNP; the lineage name of
  BOV is now "*M. bovis*, *M. caprae*, *M. orygis*".
* Consistency warnings for RD7 against lineages 5 and 6, and for lineage SNPs in *M. canettii*.
* Table: 2 new columns at the end, `RD7` and `RD12`. PDF report: RD7 and RD12 in the regions of difference table.
* SIT and SITVIT2 family: new `spoligotyper-download-sit` command, which downloads the 9,656 SITVIT2 patterns
  (3,850 SITs) published under GPL-3.0 with SpolLineages, from GitHub or its Zenodo mirror, checked with their
  checksum. New `SIT`, `SITVIT2family` and `ClosestSIT` columns, `--sit-db` option, and SIT in the PDF report (with
  the checksum of the database), the JSON and MultiQC. See [Installation](Installation#sit-database).
* [Validation](Validation) extended to 31 genomes: *M. caprae*, *M. orygis*, *M. microti*, *M. pinnipedii*,
  *M. mungi*, and more *M. africanum* and *M. canettii*.

## 0.4.2 (2026-09-24)
* PDF report: samples with the same spoligotype are grouped in the summary (largest groups first, alternate groups
  shaded, failed samples last), and each sample has its own page, in the same order, so that its tables are never
  split between pages. The review box is never separated from its heading.

## 0.4.1 (2026-09-23)
**Fixed**
* Seal failed ("Could not create the Java Virtual Machine") when a path contained "xmx" or "xms", which BBTools 40
  reads as a Java memory setting: in an input file name, or, occasionally, in the random name of the temporary folder.
* Gzipped files without the `.gz` extension, and files without a sequence extension (e.g. Galaxy's `.dat` files),
  were misread by Seal: they are now given to Seal with the extension matching their content.

## 0.4.0 (2026-09-23)
**New**
* Species check from the regions of difference RD9, RD4 and RD1: *M. tuberculosis*, *M. africanum*, *M. bovis*,
  BCG, other animal-adapted lineages, or "MTBC not detected". See [Species and lineage](Species-and-lineage).
* Lineage from the 62-SNP barcode of Coll *et al.* (2014): lineages 1 to 7 and their sublineages, *M. bovis*, with
  the lineage name and its typical spoligotype families.
* Contamination: estimated fraction of MTBC reads, with a warning below 60%.
* Mixed samples: flagged from the lineage SNPs (both alleles), incompatible lineages, partial regions of difference
  and weak spacers.
* Consistency warnings between the spoligotype, the regions of difference and the lineage.
* Closest known patterns (up to 3 spacers different) for patterns not in the database.
* JSON output with all the results and run information, and a MultiQC custom content file.
* `-j`/`--jobs`: type several samples at the same time in batch mode. `--no-species` to skip the new checks.
* PDF report: species, lineage, regions of difference and lineage SNPs for each sample; checksums of the new
  reference data; updated method and references.
* [Validation](Validation) on 16 reference genomes and 8 read sets (`validation/run_validation.sh`), and a comparison
  with SpoTyping.
* nf-core module and Galaxy tool, in `integrations/`.
* `scripts/make_reference_data.py` rebuilds the species and lineage reference data from public genomes.

**Changed**
* The table has 8 new columns at the end: `Species`, `Lineage`, `LineageName`, `RD9`, `RD4`, `RD1`,
  `MTBCFraction` and `Closest`.
* A second pass over the reads (a few seconds) for the lineage SNPs.
* The low-depth warning uses the estimated MTBC depth (total depth times the fraction of MTBC reads), shown to one
  decimal.
* Fasta files with more than 1,000 sequences and 13 Mb are typed as reads (reads in fasta format), not as an assembly.
* `--no-pdf` no longer drops the MD5 checksums from the JSON; only `--no-md5` does.
* Symbolic links to folders are followed in batch mode.

**Fixed**
* Input files, output folders or an installation path with a space or a comma made Seal fail. They are now passed
  to Seal through links with safe names.
* A single fastq file of interleaved pairs was processed by Seal as paired-end reads, doubling the counts and the
  estimated MTBC fraction.
* Seal errors are reported with their cause (e.g. "truncated or corrupt input") instead of a generic Java message.
* A sample failing with an unexpected error no longer stops a batch.
* Nanopore data: type the assembly rather than the reads, which gave inaccurate spoligotypes in our tests (see the
  [FAQ](FAQ#can-i-type-nanopore-reads)).

## 0.3.0 (2026-09-23)
**New**
* Batch mode: `-i FOLDER` types every sample of a folder and its subfolders. fastq and fasta files are detected,
  R1/R2 files are paired, and ambiguous sample names are reported rather than guessed. A sample that fails is
  reported and the others are still typed. See [Usage](Usage#a-folder-of-samples-batch-mode).
* PDF report for quality assurance: summary with the spoligotype patterns, one section per sample with the reads
  supporting each spacer, input files with size, date and MD5 checksum, and the run information (operator, user,
  computer, time zone, command, parameters, versions of spoligotyper, Python, BBTools and Java, checksums of the
  reference data), with a review and signature box. New options `--no-pdf`, `--no-md5` and `--operator`.
  See [Output files](Output-files#pdf-report-spoligotyping_reportpdf-or-sample_spoligotypingpdf).
* The table has 6 new columns after the original 6: `FileType`, `Reads`, `Depth`, `MinCount`, `Status` and
  `Warnings`.
* Warning when the estimated sequencing depth is below 20x.
* Logo.
* Zenodo DOI: https://doi.org/10.5281/zenodo.22926160 (all versions), in the README, `CITATION.cff` and the PDF report.

**Changed**
* ReportLab is now required (installed automatically with conda and pip).
* Seal errors start with their cause instead of the full command.

## 0.2.0 (2026-09-23)
Package revamp: installable with conda (bioconda) or pip, with a `spoligotyper` command.

**Fixes**
* Fixed: version 0.1 no longer started with recent versions of setuptools (`No module named 'pkg_resources'`).
* Fixed: the default `--min-count` was 4, although the help and documentation said 5. It is now 5 for fastq files.
* Fixed: fasta files needed `--min-count 1`, or spacers were missed. The minimum count is now 1 for fasta files
  automatically, and spoligotyper warns if a higher value is used with a fasta file.
* Fixed: Seal errors were hidden, and led to a confusing "file not found" error. They are now reported with Seal's
  own message.
* Fixed: Seal sized its memory from the computer's free memory (196 GB on a 512 GB server), which failed on shared
  computers and clusters. It now uses 1 GB (`--memory` to change it).
* Fixed: sample names were cut at the last `_` of any file name containing `_R1` (`Iso_R10.fasta` gave `Iso`,
  `S1_R1_001.fastq.gz` gave `S1_R1`). Only read suffixes of fastq files are removed now (`_R1`, `_R1_001`, `_1`).
* Fixed: temporary Seal files were written in the output folder, where two runs with the same sample name could
  overwrite each other.

**New**
* `spoligotyper` command, `python -m spoligotyper`, and a Python API (`spoligotyper.pipeline.spoligotype`).
* `-s`/`--sample` to set the sample name, `--db` to use another spoligotype database, `--memory`.
* The file type (fasta or fastq) is detected from the content, not from the extension.
* Checks: missing or invalid input files, `-r1` and `-r2` being the same file, fasta files given with `-r2`.
* Warnings when no spacer is found, and when spacers called absent were seen in a few reads.
* The report goes to standard output and messages to standard error.
* Documentation moved to this wiki, with a [tutorial](Tutorial) on public data, and an example script
  (`examples/run_example.sh`).
* Tests (unit tests, and end-to-end tests with Seal), continuous integration, and a bioconda recipe.

**Changed**
* `-v` now means `--verbose`; use `--version` for the version.
* `spoligotyper.py` is replaced by the `spoligotyper` command.

## 0.1
First version: `python spoligotyper.py`.
