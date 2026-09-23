#!/usr/bin/env bash
# Spoligotype three public samples and compare the results with expected_results.tsv:
#   AF2122_97  M. bovis AF2122/97, Illumina single-end reads (ENA ERR1744454, 210 MB)  -> SB0140
#   NC_002945  M. bovis AF2122/97 reference genome assembly (NCBI NC_002945.4)         -> SB0140
#   H37Rv      M. tuberculosis H37Rv reference genome assembly (NCBI NC_000962.3)      -> not an SB pattern
# Downloads go to ./data/ (kept for later runs) and results to ./results/.
# Usage: bash run_example.sh [threads]
set -euo pipefail

cd "$(dirname "$0")"
threads=${1:-4}
command -v spoligotyper >/dev/null || { echo 'spoligotyper is not installed: see the Installation page of the wiki' >&2; exit 1; }

download() {  # url output
    [ -s "$2" ] && return
    echo "Downloading $2" >&2
    curl -sSfL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
}

efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id='
mkdir -p data results
download https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR174/004/ERR1744454/ERR1744454.fastq.gz data/AF2122_97.fastq.gz
download "${efetch}NC_002945.4" data/NC_002945.fasta
download "${efetch}NC_000962.3" data/H37Rv.fasta

spoligotyper -r1 data/AF2122_97.fastq.gz -o results/ -t "$threads"
spoligotyper -r1 data/NC_002945.fasta -o results/ -t "$threads"
spoligotyper -r1 data/H37Rv.fasta -o results/ -t "$threads"

# One table for all the samples: the header once, then one line per sample
awk 'FNR == 1 && NR > 1 {next} 1' results/*_spoligotyping.txt > results/all_samples.tsv

# SpacerCount (column 2) depends on the BBTools version; the spoligotype does not
if diff <(cut -f 1,3- expected_results.tsv) <(cut -f 1,3- results/all_samples.tsv); then
    echo "OK: results match expected_results.tsv (results/all_samples.tsv)"
else
    echo "Results differ from expected_results.tsv" >&2
    exit 1
fi
