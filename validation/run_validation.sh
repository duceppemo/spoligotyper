#!/usr/bin/env bash
# Validate spoligotyper on public reference genomes (genomes.tsv), public reads, and simulated reads:
# pure, mixed (two strains) and contaminated (MTBC + M. marinum). Writes results/validation.md.
# Requires spoligotyper, BBTools (seal.sh, randomreads.sh) and curl. About 600 MB of downloads, kept in data/.
# Usage: bash run_validation.sh [threads]
set -euo pipefail

cd "$(dirname "$0")"
threads=${1:-8}
efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id='
mkdir -p data/genomes data/reads results

download() {  # url output
    [ -s "$2" ] && return
    echo "Downloading $2" >&2
    curl -sSfL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
}

# Reference genomes, named after the sample column
grep -v '^#' genomes.tsv | tail -n +2 | while IFS=$'\t' read -r sample accession _; do
    download "${efetch}${accession}" "data/genomes/${sample}.fasta"
    sleep 0.4  # NCBI: at most 3 requests per second
done

# Public reads: M. bovis AF2122/97, Illumina single-end
download https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR174/004/ERR1744454/ERR1744454.fastq.gz data/reads/ERR1744454.fastq.gz

# Simulated 150 bp reads with sequencing errors, fixed seeds. G = genome size / read length, per 1x of depth.
simulate() {  # genome depth seed output [paired]
    local out=$4
    [ -s "$out" ] && return
    local reads=$(( $(grep -v '>' "data/genomes/$1.fasta" | tr -d '\n' | wc -c) * $2 / 150 ))
    if [ "${5:-}" = paired ]; then
        randomreads.sh -Xmx2g ref="data/genomes/$1.fasta" out="$out" out2="${out/_R1/_R2}" length=150 path=data/sim \
            reads=$(( reads / 2 )) paired=t mininsert=250 maxinsert=450 seed="$3" adderrors=t ow=t 2>/dev/null
    else
        randomreads.sh -Xmx2g ref="data/genomes/$1.fasta" out="$out" length=150 path=data/sim reads="$reads" seed="$3" \
            adderrors=t ow=t 2>/dev/null
    fi
}
mkdir -p data/sim
simulate H37Rv 30 1 data/reads/sim_H37Rv_30x.fastq.gz
simulate H37Rv 30 2 data/reads/sim_H37Rv_30x_PE_R1.fastq.gz paired
simulate H37Rv 10 3 data/reads/sim_H37Rv_10x.fastq.gz
simulate CCDC5079 30 4 data/reads/sim_Beijing_30x.fastq.gz
simulate AF2122_97 30 5 data/reads/sim_AF2122_97_30x.fastq.gz
simulate H37Rv 21 6 data/sim/H37Rv_21x.fastq.gz
simulate AF2122_97 9 7 data/sim/AF2122_97_9x.fastq.gz
simulate H37Rv 15 8 data/sim/H37Rv_15x.fastq.gz
simulate M_marinum 15 9 data/sim/M_marinum_15x.fastq.gz
[ -s data/reads/sim_mixed_H37Rv70_AF2122_30.fastq.gz ] ||
    cat data/sim/H37Rv_21x.fastq.gz data/sim/AF2122_97_9x.fastq.gz > data/reads/sim_mixed_H37Rv70_AF2122_30.fastq.gz
[ -s data/reads/sim_contaminated_H37Rv15x_marinum15x.fastq.gz ] ||
    cat data/sim/H37Rv_15x.fastq.gz data/sim/M_marinum_15x.fastq.gz \
        > data/reads/sim_contaminated_H37Rv15x_marinum15x.fastq.gz

spoligotyper -i data/genomes -o results/genomes -t "$threads" -j 4 --operator validation || true
spoligotyper -i data/reads -o results/reads -t "$threads" -j 2 --operator validation || true
python3 check_results.py > results/validation.md
cat results/validation.md
