# FAQ

### My sample is `Spoligo not found`
The pattern is not in the [Mbovis.org](https://www.mbovis.org/) database. The `Closest` column lists the closest SB
numbers (up to 3 spacers different). Either:
* **It is a human-adapted lineage** (*M. tuberculosis*): SB numbers only exist for the RD9-deleted lineages
  (*M. bovis*, *M. caprae*, *M. pinnipedii*, *M. microti*, *M. africanum*, ...). The `Lineage` column gives the
  lineage and its typical spoligotype families, and the `SIT` and `SITVIT2family` columns its shared international
  type (SIT), once the [SIT database](Installation#sit-database) is downloaded (H37Rv, 777777477760771, is SIT451).
* **It is a new pattern**: new *M. bovis* patterns can be submitted to Mbovis.org to get an SB number.
* **A spacer is miscalled**: look at the `SpacerCount` column for counts close to `--min-count`, and see
  [How it works](How-it-works#minimum-count).

### A sample with no spacer is SB2277
SB2277 is the pattern with no spacer at all. A sample that is not from the *M. tuberculosis* complex, or a file with
no MTBC reads, also has no spacer, and so gets SB2277. spoligotyper flags these samples with a "no spacer found"
warning, and the `Species` column says "MTBC not detected" when there is no MTBC DNA. *M. canettii* has no standard
spacer either: it is reported with RD9 intact and no lineage.

### How reliable are the species and the lineage?
They were checked on 16 reference genomes and 8 read sets of known species and lineage: see [Validation](Validation).
The lineage comes from the SNP barcode of Coll *et al.* (2014), the reference method for SNP-based lineage typing.
For drug resistance and a finer lineage, use a dedicated tool such as [TB-Profiler](https://github.com/jodyphelan/TBProfiler).

### What does "only about X% of the reads appear to be from the M. tuberculosis complex" mean?
The MTBC-specific control regions have fewer reads than the sequencing depth predicts: part of the reads come from
something else (contamination, host DNA, another organism). The spoligotype and species are still called from the
MTBC reads. See [Species and lineage](Species-and-lineage#contamination-fraction-of-mtbc-reads).

### Where do the SIT numbers come from?
From the SITVIT2 database of the Institut Pasteur de Guadeloupe, through the 9,656 SITVIT2 patterns published under
GPL-3.0 with SpolLineages (3,850 SITs, up to SIT3862, as of 2022). `spoligotyper-download-sit` downloads them: see
[Installation](Installation#sit-database). SITs created since 2022 are missing: a pattern without SIT is reported as
`Orphan` (known to SITVIT2, seen once) or `Spoligo not found`, with the closest SITs. SITVIT2 families are given for
orphan patterns too.

### Can I use another database?
Yes, with `--db my_database.txt`. The file has one pattern per line, with 3 columns separated by spaces or tabs:
the octal code, the name, and the binary pattern. Lines starting with `#` are ignored. For example, with SIT numbers:
```
# octal          name    binary
777777477760771  SIT451  1111111111111111111001111111111100001111111
```
spoligotyper checks that each octal code matches its binary pattern.

### How recent is the included database?
It is a snapshot of the Mbovis.org database with 1,976 SB patterns. To use a newer version, download it from
Mbovis.org, format it as above and use `--db`.

### Can I type nanopore reads?
Type the **assembly**, not the reads. In our hands, spoligotypes called directly from nanopore reads were often
wrong, while spoligotypes called from assemblies of the same long reads were almost always right.

Spacers are 25 bp long and spoligotyper allows 1 mismatch per spacer: individual nanopore reads carry enough errors
(substitutions and, above all, small insertions and deletions) that many reads covering a present spacer are missed,
and present spacers can fall below `--min-count`. The consensus sequence of an assembly corrects these errors. So:
1. assemble the long reads (e.g. [Flye](https://github.com/mikolmogorov/Flye) or
   [Autocycler](https://github.com/rrwick/Autocycler)), ideally with polishing (e.g.
   [Medaka](https://github.com/nanoporetech/medaka));
2. type the assembly: `spoligotyper -r1 assembly.fasta -o results/`.

The direct repeat locus is repetitive, but long reads usually span it entirely, so it assembles well.
If you type nanopore reads anyway, treat the result as provisional: check the `SpacerCount` column for present
spacers with low counts, and confirm with the assembly.

### Does spoligotyper work on *M. tuberculosis*?
Yes: spacer detection and the binary, octal and hexadecimal codes work for the whole complex. Only the SB number
is specific to the animal-adapted lineages (see above).

### Why are there 43 spacers?
The 43 spacers are those of the standard spoligotyping membrane
([Kamerbeek *et al.* 1997](https://doi.org/10.1128/jcm.35.4.907-914.1997)), so that *in silico* results can be
compared with laboratory spoligotyping. Other spacers exist in some strains but are not part of the standard
pattern.

### How do I type many samples?
Put them in a folder and use `-i`: see [Usage](Usage#a-folder-of-samples-batch-mode). fastq and fasta files are
detected, and R1/R2 files are paired automatically. You get one table and one PDF report for the whole run.

### Can the PDF report be used for accredited (ISO 17025) work?
It is designed for it: it records the operator, date and time with time zone, computer, exact command, parameters,
versions of all the software, checksums of the input files and of the reference data, and the evidence behind each
call, and it has a review and signature box. Validating the method for your scope remains your laboratory's
responsibility; the [tutorial](Tutorial) data, with known spoligotypes, can be part of it.

### Should I type reads or the assembly?
Reads, when you have them: the DR locus is repetitive and can be broken or collapsed in short-read assemblies.
See [How it works](How-it-works#limitations).

### Why are the times in the PDF report in UTC?
spoligotyper writes every date and time of the report (start, end, file dates, footer) in the time zone of the
computer that runs it, and shows that time zone (e.g. `EDT`, `UTC`). Containers usually run in UTC: with Docker, give
the container the time zone of the host with `-v /etc/localtime:/etc/localtime:ro` (or `-e TZ=America/Toronto` if the
image has time zone data); in Nextflow, add it to `docker.runOptions`. Apptainer/Singularity containers normally use
the time zone of the host already.

### How do I cite spoligotyper?
See [Home](Home#citing).
