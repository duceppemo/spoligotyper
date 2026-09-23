# Installation

## Requirements
* Linux or macOS
* Python 3.10 or later (no other Python package needed)
* [BBTools](https://sourceforge.net/projects/bbmap/) 38 or later, for `seal.sh`, and Java (installed with BBTools by
  conda)

## With conda (recommended)
```
conda create -n spoligotyper -c conda-forge -c bioconda spoligotyper
conda activate spoligotyper
spoligotyper --version
```
This installs BBTools and Java too.

## With pip
```
pip install spoligotyper
spoligotyper --version
```
BBTools is not available from PyPI and must be installed separately, for example with
`conda install -c bioconda bbmap`, or from the [BBTools downloads](https://sourceforge.net/projects/bbmap/).
`spoligotyper` finds `seal.sh` in your `PATH`, or next to the Python interpreter it runs with.

## From source
```
git clone https://github.com/duceppemo/spoligotyper
cd spoligotyper
conda env create -f environment.yml
conda activate spoligotyper
pip install --no-deps .
spoligotyper --version
```

## Without installing
From the cloned folder, as long as `seal.sh` is in your `PATH`:
```
python3 -m spoligotyper -h
```

## Updating from 0.1
Version 0.2.0 is an installable package with a `spoligotyper` command; the `spoligotyper.py` script is gone.
Replace `python spoligotyper.py ...` by `spoligotyper ...` (or `python -m spoligotyper ...` from the cloned folder).
The options are the same, except:
* `-v` now means `--verbose`; use `--version` for the version.
* `--min-count` no longer needs to be set for fasta files: it is 1 for fasta files and 5 for fastq files
  (0.1 used 4, although its help said 5).

The report has the same name and columns. See the [Changelog](Changelog) for all the changes.
