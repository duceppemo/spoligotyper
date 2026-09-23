# Changelog

## 0.3.0 (2026-09-24)
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
