<p align="center"><img src="https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/logo.svg" width="520" alt="spoligotyper logo"></p>

# spoligotyper

**spoligotyper** performs *in silico* spoligotyping of *Mycobacterium tuberculosis* complex (MTBC) samples from
sequencing reads (fastq) or genome assemblies (fasta). It finds the 43 spacers of the direct repeat (DR) locus with
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools, and reports the spoligotype as binary, octal and
hexadecimal codes, and as an SB number from the [Mbovis.org](https://www.mbovis.org/) database.

It types one sample or a whole folder of samples (fastq and fasta files are detected and R1/R2 files paired
automatically), and writes a table and a [PDF report](Output-files#pdf-report-spoligotyping_reportpdf-or-sample_spoligotypingpdf)
with everything needed for quality assurance: the evidence for each call, input files and checksums, software
versions, parameters, operator and date.

It works for all members of the complex (*M. tuberculosis*, *M. bovis*, *M. caprae*, *M. pinnipedii*,
*M. microti*, *M. africanum*, *M. canettii*, ...). SB numbers are only defined for the animal-adapted lineages
(RD9-deleted: *M. bovis*, *M. caprae*, ...); for human-adapted lineages, use the octal code (see the [FAQ](FAQ)).

![Summary page of the PDF report](https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/report_summary.png)

## Getting started
* [Installation](Installation)
* [Tutorial](Tutorial): three public genomes, step by step
* [Usage](Usage): input files, options and examples

## Understanding the results
* [Output files](Output-files): the table, column by column, and the PDF report
* [How it works](How-it-works): spoligotyping, spacer detection, and the choice of the minimum count

## Reference
* [Troubleshooting](Troubleshooting)
* [FAQ](FAQ)
* [Changelog](Changelog)
* [Contributing](Contributing)

## Citing
Please cite spoligotyper: Duceppe M-O. spoligotyper: in silico spoligotyping of Mycobacterium tuberculosis complex
genomes. Zenodo. https://doi.org/10.5281/zenodo.22926160 (all versions; each release also has its own DOI on Zenodo). Please also cite
BBTools: Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/
