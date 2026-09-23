<h1 align="center">spoligotyper</h1>

<p align="center">
  <a href="https://github.com/duceppemo/spoligotyper/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/spoligotyper/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://codecov.io/gh/duceppemo/spoligotyper"><img src="https://codecov.io/gh/duceppemo/spoligotyper/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/spoligotyper/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/spoligotyper" alt="Release"></a>
  <a href="https://pypi.org/project/spoligotyper/"><img src="https://img.shields.io/pypi/v/spoligotyper" alt="PyPI"></a>
  <a href="https://bioconda.github.io/recipes/spoligotyper/README.html"><img src="https://img.shields.io/conda/vn/bioconda/spoligotyper?label=bioconda" alt="Bioconda"></a>
  <img src="https://img.shields.io/badge/python-3.10%E2%80%933.14-blue" alt="Python 3.10–3.14">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/spoligotyper" alt="License"></a>
  <a href="https://github.com/duceppemo/spoligotyper/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
</p>

*In silico* spoligotyping of *Mycobacterium tuberculosis* complex (MTBC) samples, from sequencing reads (fastq) or
genome assemblies (fasta). spoligotyper finds the 43 spacers of the direct repeat (DR) locus with
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools and reports the spoligotype as binary, octal and
hexadecimal codes, and as an SB number from the [Mbovis.org](https://www.mbovis.org/) database.

```
Sample     SpacerCount           Binary                                       Octal            Hexadecimal        Spoligotype
AF2122_97  56:47:0:58:63:0:...   1101101000001110111111111111111111111100000  664073777777600  6D-03-5F-7F-FF-60  SB0140
```

## Features
* **Reads or assemblies**: single-end or paired-end fastq, or fasta, gzipped or not. Illumina and nanopore reads.
* **Fast**: a few seconds per sample, with 1 GB of memory.
* **All the standard codes**: binary, octal, hexadecimal, and SB number for *M. bovis* and other animal-adapted lineages.
* **Transparent**: the number of reads supporting each spacer is reported, and borderline calls are flagged.

## Installation
```
conda install -c conda-forge -c bioconda spoligotyper
```
Or with pip, if BBTools is already installed (`conda install -c bioconda bbmap`): `pip install spoligotyper`.
See [Installation](https://github.com/duceppemo/spoligotyper/wiki/Installation) for other options.

## Quick start
```
spoligotyper -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/   # Paired-end reads
spoligotyper -r1 sample.fastq.gz -o results/                             # Single-end or nanopore reads
spoligotyper -r1 assembly.fasta -o results/                              # Assembly
```
The report is printed and saved as `results/<sample>_spoligotyping.txt`. New to the tool? The
**[tutorial](https://github.com/duceppemo/spoligotyper/wiki/Tutorial)** types three public genomes in a few minutes.

## Documentation
The **[wiki](https://github.com/duceppemo/spoligotyper/wiki)** covers
[usage and options](https://github.com/duceppemo/spoligotyper/wiki/Usage),
[output files](https://github.com/duceppemo/spoligotyper/wiki/Output-files),
[how it works](https://github.com/duceppemo/spoligotyper/wiki/How-it-works),
[troubleshooting](https://github.com/duceppemo/spoligotyper/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/spoligotyper/wiki/FAQ).

## Citing
If you use spoligotyper, please cite this repository (GitHub's **"Cite this repository"** button gives APA and
BibTeX formats), and BBTools: Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/

## Contributing
Bug reports, questions and pull requests are welcome: see [CONTRIBUTING.md](CONTRIBUTING.md).

## Author
Marc-Olivier Duceppe, Canadian Food Inspection Agency (CFIA): marc-olivier.duceppe@inspection.gc.ca

## License
[MIT](LICENSE)
