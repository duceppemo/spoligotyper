# Installation

## Requirements
* Linux or macOS
* Python 3.10 or later, with [ReportLab](https://www.reportlab.com/opensource/) 4 or later (for the PDF report)
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

## SIT database
SB numbers (Mbovis.org) are included with spoligotyper. Shared international types (SIT) and SITVIT2 families come
from the SITVIT2 database, which is not openly licensed and cannot be included. Instead, spoligotyper uses the 9,656
SITVIT2 patterns (3,850 with a SIT) published under GPL-3.0 with [SpolLineages](https://github.com/dcouvin/SpolLineages)
(Couvin *et al.* 2020). Download them once:
```
spoligotyper-download-sit
```
The list is downloaded from the SpolLineages repository (a fixed version, checked with its SHA-256 checksum), or from
its Zenodo mirror if GitHub cannot be reached, and saved in `~/.cache/spoligotyper/` (or in the folder of the
`SPOLIGOTYPER_DATA` environment variable). spoligotyper then fills the `SIT`, `SITVIT2family` and `ClosestSIT`
columns automatically.

* **Computers without internet access** (e.g. cluster nodes): run `spoligotyper-download-sit -o /shared/folder` on a
  computer with access, then `spoligotyper --sit-db /shared/folder/sit_database.tsv ...`, or set
  `SPOLIGOTYPER_DATA=/shared/folder`. Alternatively, copy `Spoligo_list.csv` and run
  `spoligotyper-download-sit --source Spoligo_list.csv`.
* The list dates from 2022: SITs created in SITVIT2 since then are missing (the `ClosestSIT` column lists the closest
  known ones).
* Please cite Couvin *et al.* 2020 (https://doi.org/10.1093/database/baaa108) and Couvin *et al.* 2019
  (https://doi.org/10.1016/j.meegid.2018.12.030) when you report SITs.

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

## Updating from 0.3
Version 0.4.0 adds the species check and the lineage, a JSON file and a MultiQC file. The table has 8 new columns
after the 12 of version 0.3. Use `--no-species` for the previous behaviour (one pass over the reads).

## Updating from 0.2
Version 0.3.0 needs ReportLab: `conda install -c conda-forge reportlab` or `pip install reportlab` (installed
automatically with conda or pip). The table has 6 new columns after the original 6, and a PDF report is written
next to it (`--no-pdf` to skip it). See the [Changelog](Changelog).

## Updating from 0.1
Version 0.2.0 is an installable package with a `spoligotyper` command; the `spoligotyper.py` script is gone.
Replace `python spoligotyper.py ...` by `spoligotyper ...` (or `python -m spoligotyper ...` from the cloned folder).
The options are the same, except:
* `-v` now means `--verbose`; use `--version` for the version.
* `--min-count` no longer needs to be set for fasta files: it is 1 for fasta files and 5 for fastq files
  (0.1 used 4, although its help said 5).

The report has the same name and columns. See the [Changelog](Changelog) for all the changes.
