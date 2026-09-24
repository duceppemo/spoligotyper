<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/logo_dark.svg">
    <img src="https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/logo.svg" width="520" alt="spoligotyper: in silico spoligotyping of the M. tuberculosis complex">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/duceppemo/spoligotyper/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/spoligotyper/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://codecov.io/gh/duceppemo/spoligotyper"><img src="https://codecov.io/gh/duceppemo/spoligotyper/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/spoligotyper/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/spoligotyper?cacheSeconds=3600" alt="Release"></a>
  <a href="https://pypi.org/project/spoligotyper/"><img src="https://img.shields.io/pypi/v/spoligotyper?cacheSeconds=3600" alt="PyPI"></a>
  <a href="https://bioconda.github.io/recipes/spoligotyper/README.html"><img src="https://img.shields.io/conda/vn/bioconda/spoligotyper?label=bioconda&cacheSeconds=3600" alt="Bioconda"></a>
  <img src="https://img.shields.io/badge/python-3.10%E2%80%933.14-blue" alt="Python 3.10–3.14">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/spoligotyper" alt="License"></a>
  <a href="https://github.com/duceppemo/spoligotyper/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.22926160"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.22926160.svg" alt="DOI"></a>
</p>

*In silico* spoligotyping of *Mycobacterium tuberculosis* complex (MTBC) samples, from sequencing reads (fastq) or
genome assemblies (fasta). spoligotyper finds the 43 spacers of the direct repeat (DR) locus with
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools and reports the spoligotype as binary, octal and
hexadecimal codes, and as an SB number from the [Mbovis.org](https://www.mbovis.org/) database. It also identifies
the species and the lineage, and flags contaminated or mixed samples.

<p align="center">
  <img src="https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/report_summary.png" width="49%" alt="Summary page of the PDF report, with the spoligotype pattern of three samples">
  <img src="https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/report_sample.png" width="49%" alt="Sample page of the PDF report, with the reads supporting each spacer">
  <br>
  <sub>PDF report of the <a href="https://github.com/duceppemo/spoligotyper/wiki/Tutorial">tutorial</a>: summary, and the evidence for each sample.</sub>
</p>

## Features
* **Reads or assemblies**: single-end or paired-end fastq, or fasta, gzipped or not. For nanopore data, type the
  assembly ([why](https://github.com/duceppemo/spoligotyper/wiki/FAQ#can-i-type-nanopore-reads)).
* **Batch mode**: point it at a folder; fastq and fasta files are detected and R1/R2 files paired automatically.
* **All the standard codes**: binary, octal, hexadecimal, and SB number for *M. bovis* and other animal-adapted lineages.
* **Species and lineage**: *M. tuberculosis*, *M. africanum*, *M. bovis*, BCG, ... from regions of difference, and
  the lineage (1 to 7 and sublineages) from a 62-SNP barcode.
* **Quality checks**: fraction of MTBC reads (contamination), mixed samples, consistency between the spoligotype,
  species and lineage, and the closest known patterns for new spoligotypes.
* **PDF report for QA**: results, reads supporting each spacer, input files with checksums, software versions,
  parameters, operator, date, and a review box. Plus a table, JSON and a MultiQC section for pipelines.
* **Validated** on 16 reference genomes and 8 read sets of known species, lineage and spoligotype
  ([Validation](https://github.com/duceppemo/spoligotyper/wiki/Validation)).
* **Workflow ready**: nf-core module and Galaxy tool in [`integrations/`](integrations/).
* **Transparent**: borderline calls, low depth and failed samples are flagged, never hidden.
* **Fast**: a few seconds per sample, with 1 GB of memory.

## Installation
```
conda install -c conda-forge -c bioconda spoligotyper
```
Or with pip, if BBTools is already installed (`conda install -c bioconda bbmap`): `pip install spoligotyper`.
See [Installation](https://github.com/duceppemo/spoligotyper/wiki/Installation) for other options.

## Quick start
```
spoligotyper -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/   # Paired-end reads
spoligotyper -r1 sample.fastq.gz -o results/                             # Single-end reads
spoligotyper -r1 assembly.fasta -o results/                              # Assembly (also for nanopore data)
spoligotyper -i folder/ -o results/                                      # All the samples in a folder
```
The results are printed and saved in `results/` as a table and a PDF report. New to the tool? The
**[tutorial](https://github.com/duceppemo/spoligotyper/wiki/Tutorial)** types three public genomes in a few minutes.

## Documentation
The **[wiki](https://github.com/duceppemo/spoligotyper/wiki)** covers
[usage and options](https://github.com/duceppemo/spoligotyper/wiki/Usage),
[output files](https://github.com/duceppemo/spoligotyper/wiki/Output-files),
[how it works](https://github.com/duceppemo/spoligotyper/wiki/How-it-works),
[species and lineage](https://github.com/duceppemo/spoligotyper/wiki/Species-and-lineage),
[validation](https://github.com/duceppemo/spoligotyper/wiki/Validation),
[troubleshooting](https://github.com/duceppemo/spoligotyper/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/spoligotyper/wiki/FAQ).

## Citing
If you use spoligotyper, please cite it and BBTools, which finds the spacers:

> Duceppe M-O. spoligotyper: in silico spoligotyping of Mycobacterium tuberculosis complex genomes. Zenodo.
> https://doi.org/10.5281/zenodo.22926160

This DOI always points to the latest version; each release also has its own DOI, listed on
[Zenodo](https://doi.org/10.5281/zenodo.22926160). GitHub's **"Cite this repository"** button gives the same citation in APA and BibTeX
formats.

> Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/

For the lineage, please also cite the SNP barcode:

> Coll F *et al.* A robust SNP barcode for typing *Mycobacterium tuberculosis* complex strains. *Nat Commun* 5, 4812
> (2014). https://doi.org/10.1038/ncomms5812

## Contributing
Bug reports, questions and pull requests are welcome: see [CONTRIBUTING.md](CONTRIBUTING.md).

## Author
Marc-Olivier Duceppe, Canadian Food Inspection Agency (CFIA): marc-olivier.duceppe@inspection.gc.ca

## License
[MIT](LICENSE)
