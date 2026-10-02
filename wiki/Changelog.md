# Changelog

## Unreleased
* **Lineage 1 sublineages** after [Netikul *et al.* 2022](https://doi.org/10.1038/s41598-022-05524-0): L1.1 to
  L1.3 and their sublineages down to the fourth level (32 in all), from the 1,835 sublineage-specific SNPs of the
  paper, in the same pass as the other SNPs. New `L1Sublineage` column (at the end), MultiQC column, JSON `l1`, and a
  line in the PDF report with the typical spoligotype families. The names are those of the revised lineage 1
  nomenclature, kept apart from those of Coll *et al.* in the `Lineage` column (Coll 1.2.1 = L1.2, Coll 1.2.2 = L1.3,
  Coll 1.1.2 = L1.1.2.2). A sublineage is called when at least 2 of its SNPs and half of those covered carry the
  derived allele, and its parent group is called; mixed when 2 of its SNPs and 10% of those covered have both
  alleles. Warning when lineage 1 sublineage SNPs are found in a sample without lineage 1 SNPs of Coll *et al.*
* A region partially deleted in one strain no longer makes the RD profile "unusual": when the profile is unknown, a
  partially deleted region counts as present if that gives a known profile (e.g. a lineage 1 strain lacking 9 kb of
  RD7 stays *M. tuberculosis*; the partial deletion is still reported).
* **PDF report**, reorganized: two pages per sample, (1) the spoligotype, input and reads per spacer, with their
  warnings, (2) the species and lineage, with tables of the livestock lineage and lineage 1 sublineage groups (SNPs
  covered, with the derived allele, with both alleles, reads, result) and their warnings; a page too long (long input
  paths, many warnings) is scaled down to fit. Summary and warnings tables: every other row shaded, groups of samples
  with the same spoligotype between blue lines, a header row for the warnings, and tables over several pages titled
  "(page 1 of n)" with their header repeated. Short RD results (the coordinates of partial deletions are in the
  warnings). Checked on the 84 samples of the validation (176 pages).
* Validation (84 of 84 checks): one public run of each of the 25 terminal lineage 1 sublineages, from the genomes of
  the paper.
* **Dassie bacillus**: identified from its own partial deletion of RD1, RD1<sup>das</sup> (Rv3874 to Rv3877), in the
  row of *M. microti* (RD1 segments 4 to 12 missing; 11 of 20 found, counted as +). Before, its row required RD1
  deleted, which its RD1<sup>das</sup> never gives: it would have been reported as "*M. microti*, *M. pinnipedii* or
  *M. mungi*". Not validated: no genome is public.
* Lineages 4 and 4.9 (both the H37Rv allele, which some *M. canettii* carry) no longer make a strain without
  standard spacers *M. tuberculosis* (4.9 did).
* **File names**: only paths of plain characters are given to `seal.sh` (which runs its command line through the
  shell); others (quotes, `$`, backquotes, `;`, `&`, parentheses...) are linked under neutral names. A file name
  could run a command. A fastq file named `.fasta` (or the reverse) is now read as its content says.
* Sample names: tabs and line breaks are rejected (`--sample`) or replaced. SIT database: blank lines are ignored.
* Warnings: "no spacer found" when no spacer is present, even with a few reads on some (with a specific text for
  MTBC samples: deletion of the direct repeat locus?); no "is an SB number" warning for the pattern without spacers.
* Data: spacers 1 and 2 in the standard orientation (no effect: both strands are searched); spoligotype families of
  lineage 4.6.2 without a double space; La1.8 includes unknown8.
* Docs: Dassie bacillus and RD1<sup>mic</sup> (Brodin *et al.* 2002, Mostowy *et al.* 2004 cited), lineage 4.9,
  installation with pip until the bioconda package is available, no Zenodo mirror of the SIT list, Validation page
  claims corrected, SpoTyping results of the comparison in `validation/spotyping.tsv`.
* The livestock lineages and the lineage 1 sublineages share one implementation (`snp_groups`). JSON: `livestock`
  and `l1` have two new keys, `main` and `mixed_within`, and only the SNPs of the groups with a derived allele (not
  all 1,835 for every sample). Mixed texts list at most 6 groups.

## 0.8.0 (2026-10-02)
* **Livestock lineages** after [Zwyer *et al.* 2021](https://doi.org/10.12688/openreseurope.14029.2): La1
  (*M. bovis*), La2 (*M. caprae*), La3 (*M. orygis*) and the La1 sublineages La1.1 to La1.8 (e.g. La1.8.1, formerly
  Eu1; BCG in La1.2), from 88 of the marker SNPs of the paper, in the same pass as the lineage SNPs. New `LaLineage`
  column (at the end of the table), MultiQC column, JSON `livestock`, and a line in the PDF report with the SNPs
  supporting each group.
* *M. caprae* and *M. orygis*, which have the same RD profile, are now told apart (La2 and La3). Warnings when the
  livestock lineage contradicts the regions of difference.
* Mixed samples: a sublineage SNP with both alleles is a mix only if the SNP of its parent lineage also has the
  lineage allele (in at least 15% of its reads). A lone lineage 2.1 SNP at 23% in an *M. bovis* strain no longer
  makes it "mixed". A livestock lineage group is mixed when at least 2 of its SNPs have both alleles, and a mix of
  sublineages of one lineage keeps its species (e.g. "M. bovis, mixed sample?").
* Spacer 3 of *M. orygis*, a variant with 2 mismatches detected by the spoligotyping membrane, is counted for spacer 3
  (`spacer_variants.fasta`), and the sample page says so: *M. orygis* now gets its Mbovis.org name (e.g. SB0422).
  The variant is only found in *M. orygis* among the validation genomes and reads (and in lineage 6, whose spacer 3
  is between the two sequences and already present).
* Validation (59 of 59 checks): one public run of each livestock lineage and sublineage, from the genomes of the
  paper, with its sublineage and SB number. The validation script downloads complete runs when they are small, and
  whole runs whose reads are sorted by position.

## 0.7.0 (2026-09-28)
* Patterns with an SB number get no closest SITs: SITVIT2 lacks many patterns of the animal-adapted lineages (only 603
  of the 1,976 Mbovis.org patterns have a SIT), for which the SB number is the reference name. The PDF report says so.
* Warning for long reads (mean length above 1,000 bp): spoligotypes from nanopore reads are often wrong; type an
  assembly of the reads.
* *M. microti* is recognised from real reads: its RD1<sup>mic</sup> deletion tolerates one error on each side (a few
  stray reads on one RD1<sup>mic</sup> segment, no read on one GC-rich RD1 segment), as in old Illumina GAII reads.
* Validation on real reads (45 of 45 checks): *M. tuberculosis* H37Rv, *M. microti*, *M. orygis*, *M. africanum*
  (lineage 6) and *M. canettii* (Illumina and nanopore), besides *M. bovis* AF2122/97. It shows what simulated reads
  cannot: an H37Rv lab stock that lost spacer 40 (SIT1647), flagged as a mixed sample for the minority that kept
  spacers 41 to 43, and a minority population of *M. orygis* with a spacer 3 variant, also flagged.

## 0.6.0 (2026-09-28)
Clearer terminology and RD reporting, after feedback from a tuberculosis expert.
* **Renamed columns**: `Spoligotype` is now `SB` and `Closest` is now `ClosestSB` in the TSV (same positions), `SB` in
  the MultiQC table, and `sb` in the JSON (`closest` entries too): they hold the name of the pattern in the Mbovis.org
  database, not the spoligotype. Scripts that read these columns need the new names.
* A pattern without SB number is reported as `Not in Mbovis.org` (`Not in database` with `--db`) instead of
  "Spoligo not found": the spoligotype is the binary, octal and hexadecimal codes, which are universal; the SB number is
  only a name in one database. Likewise, `Not in SITVIT2 list` for SITs.
* PDF report: the spoligotype codes come first (octal code in the summary and on the sample pages), and the SB number
  and SIT are labelled with their database. New "Definitions and methods" section: what the codes, names and RD
  results mean, how the RDs were determined (H37Rv coordinates, segments, thresholds), and how they compare with PCR
  assays (presence of the region's DNA, not amplicon sizes).
* Regions of difference are called from their segments: `present`, `deleted`, `partial` (some segments missing, with
  their H37Rv coordinates) or `reduced` (present at low depth: mixed sample?). Before, a median depth hid partial
  deletions, such as the part of RD1 lost by *M. microti* (RD1<sup>mic</sup>), now reported, and used to identify
  *M. microti*. In assemblies, a region is present only if all its segments are found.
  In reads, a segment is found at 5% of the control depth (GC-rich RD1 segments can drop to 10% in real Illumina
  reads), and one missing or stray segment is tolerated in the short RD9 region. Assemblies are never "reduced".
* The "is an SB number, but RD9 is present" warning is only given for SB numbers, not for names from `--db`.
* JSON: segments found, missing coordinates and H37Rv region of each RD.

## 0.5.1 (2026-09-25)
* PDF report: the date under the title now shows its time zone, like the other dates of the report. In containers,
  which usually run in UTC, see the [FAQ](FAQ#why-are-the-times-in-the-pdf-report-in-utc) to use the local time zone.

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
