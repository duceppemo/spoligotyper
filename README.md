# Spoligotyper
## Description
*In silico* spoligotyper for Mycobacterium bovis from NGS sequencing (fastq) and assembly (fasta) data. Designed to run on Linux. Tested on Ubuntu 18.

## Workflow
1. Single- or paired-end fastq files or fasta assembly are searched for the spoligo spacers using `Seal` from `BBtools`.
2. Spacer counts are converted to binary sequence (Default minimum count to be considered present is 5).
3. Binary string is converted to octal string.
4. Binary string is converted to Spoligotype.
5. Report is printed.
## Installation
1. Create and activate `spoligotyper` virtual environment:
```commandline
# Create virtual environment
conda create -n spoligotyper -y -c bioconda python=3.9 bbmap

# Activate virtual environment
conda activate spoligotyper
```
2. Clone the `spoligotyper` repository and make sure it runs smoothly:
```commandline
# Clone repo
git clone https://github.com/duceppemo/spoligotyper

# Go into cloned repo
cd spoligotyper

# Test spoligotyper
python spoligotyper -h
```
## Usage
```commandline
usage: python spoligotyper.py [-h] -r1 /path/to/sample_R1.[fastq|fasta] [-r2 /path/to/R2/fastq] [-m 5] -o /path/to/output/folder/ [-t 16] [-v]

In silico spoligotyping from WGS data for Mycobacterium bovis.

optional arguments:
  -h, --help            show this help message and exit
  -r1 /path/to/sample_R1.[fastq|fasta]
                        R1 fastq file from paired-end or single end sequencing data. Can be a fasta file. Gzipped or not. Mandatory.
  -r2 /path/to/R2/fastq
                        R2 fastq file from paired-end. Optional.
  -m 5, --min-count 5   Minimum count to consider a spacer to be present. Set to 1 if "-r1" is a fasta file. Default 5. Optional.
  -o /path/to/output/folder/, --output /path/to/output/folder/
                        Folder to hold the result files. Mandatory.
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional
  -v, --version         show program's version number and exit

```
## Examples
1. Single-end input fastq file:
```commandline
python spoligotyper.py \
    -r1 ERR1744454_1.fastq.gz \
    -o spoligotyper_output/ \
    -t 16
```
2. Paired-end fastq input file:
```commandline
python spoligotyper.py \
    -r1 ERR2747598_R1.fastq.gz \
    -r2 ERR2747598_R2.fastq.gz
    -o spoligotyper_output/ \
    -t 16
```
3. Fasta input file (genome assembly):
```commandline
# It's important to set "--min-count" to 1 for fasta files
python spoligotyper.py \
    -r1 NC_002945.4.fasta \
    --min-count 1 \
    -o spoligotyper_output/ \
    -t 16
```
