#!/usr/bin/env bash
# Validate spoligotyper on public reference genomes (genomes.tsv), public reads, and simulated reads:
# pure, mixed (two strains) and contaminated (MTBC + M. marinum). Writes results/validation.md.
# Requires spoligotyper, BBTools (seal.sh, randomreads.sh), curl and unzip. About 8 GB of downloads, kept in data/.
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

assembly() {  # accession output: all the sequences of an NCBI assembly (GCF_/GCA_)
    [ -s "$2" ] && return
    echo "Downloading $2" >&2
    local tmp
    tmp=$(mktemp -d)
    curl -sSfL --retry 3 -o "$tmp/genome.zip" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$1/download?include_annotation_type=GENOME_FASTA"
    unzip -q "$tmp/genome.zip" -d "$tmp"
    cat "$tmp"/ncbi_dataset/data/"$1"/*.fna > "$2.part" && mv "$2.part" "$2"
    rm -rf "$tmp"
}

# Reference genomes, named after the sample column
grep -v '^#' genomes.tsv | tail -n +2 | while IFS=$'\t' read -r sample accession _; do
    case "$accession" in
        GC[AF]_*) assembly "$accession" "data/genomes/${sample}.fasta" ;;
        *) download "${efetch}${accession}" "data/genomes/${sample}.fasta" ;;
    esac
    sleep 0.4  # NCBI: at most 3 requests per second
done

reads_head() {  # url output reads: the first reads of a fastq.gz file, to limit the download
    [ -s "$2" ] && return
    echo "Downloading the first $3 reads of $1" >&2
    { curl -sSfL "$1" 2>/dev/null && echo 0 > "$2.status" || echo $? > "$2.status"; } |
        { gzip -dc 2>/dev/null || true; } | head -n $(( $3 * 4 )) | gzip > "$2.part"
    # curl stops with an error (23) when head has enough reads: then check the number of reads. A complete download
    # (status 0) is a run with fewer reads
    local lines status
    lines=$(gzip -dc "$2.part" | wc -l)
    status=$(cat "$2.status")
    rm -f "$2.status"
    if [ "$lines" -ne $(( $3 * 4 )) ] && { [ "$status" -ne 0 ] || [ "$lines" -eq 0 ] || [ $(( lines % 4 )) -ne 0 ]; }; then
        rm -f "$2.part"
        echo "Download failed or incomplete: $1" >&2
        return 1
    fi
    mv "$2.part" "$2"
}

# Public reads of strains of known spoligotype, species and lineage (ENA)
ena=https://ftp.sra.ebi.ac.uk/vol1/fastq
# M. bovis AF2122/97, Illumina single-end
download $ena/ERR174/004/ERR1744454/ERR1744454.fastq.gz data/reads/ERR1744454.fastq.gz
# M. tuberculosis H37Rv, Illumina HiSeq 4000 paired-end (PRJNA634239)
for mate in 1 2; do
    download $ena/SRR120/063/SRR12006063/SRR12006063_$mate.fastq.gz data/reads/SRR12006063_$mate.fastq.gz
done
# M. microti Maus IV (the M_microti_MausIV genome), Illumina GAII paired-end (PRJEB2091), first 1.5 M pairs
for mate in 1 2; do
    reads_head $ena/ERR027/ERR027297/ERR027297_$mate.fastq.gz data/reads/ERR027297_$mate.fastq.gz 1500000
done
# M. orygis 51145 (the M_orygis_51145 genome), Illumina MiniSeq paired-end
for mate in 1 2; do
    download $ena/SRR166/049/SRR16643349/SRR16643349_$mate.fastq.gz data/reads/SRR16643349_$mate.fastq.gz
done
# M. africanum RB30001 (lineage 6, the M_africanum_RB30001 genome), Illumina HiSeq 2500 paired-end, first 1 M pairs
for mate in 1 2; do
    reads_head $ena/ERR238/008/ERR2383628/ERR2383628_$mate.fastq.gz data/reads/ERR2383628_$mate.fastq.gz 1000000
done
# M. canettii ET1291 (the M_canettii_ET1291 genome): Illumina NextSeq paired-end, and nanopore (first 30,000 reads)
for mate in 1 2; do
    download $ena/SRR186/082/SRR18636082/SRR18636082_$mate.fastq.gz data/reads/SRR18636082_$mate.fastq.gz
done
reads_head $ena/SRR230/063/SRR23035463/SRR23035463_1.fastq.gz data/reads/SRR23035463.fastq.gz 30000

runs() {  # list folder: download the runs of a list, from its columns run, first_reads and fastq
    mkdir -p "$2"
    awk -F'\t' '/^#/ {next} !header {for (i = 1; i <= NF; i++) column[$i] = i; header = 1; next}
                {print $column["run"] "\t" $column["first_reads"] "\t" $column["fastq"]}' "$1" |
    while IFS=$'\t' read -r run first fastq; do
        local mate=1
        for url in ${fastq//;/ }; do
            if [ "$first" = all ]; then
                download "$url" "$2/${run}_$mate.fastq.gz"
            else
                reads_head "$url" "$2/${run}_$mate.fastq.gz" "$first"
            fi
            mate=$(( mate + 1 ))
        done
    done
}
# Livestock lineages: one run per lineage and sublineage of Zwyer et al. 2021 (first 600,000 read pairs, or all)
runs livestock_reads.tsv data/la_reads
# Lineage 1 sublineages: one run per terminal sublineage of Netikul et al. 2022
runs l1_reads.tsv data/l1_reads

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

# SIT database (SITVIT2 patterns of SpolLineages), kept in data/
[ -s data/sit/sit_database.tsv ] || spoligotyper-download-sit -o data/sit

spoligotyper -i data/genomes -o results/genomes -t "$threads" -j 4 --operator validation --sit-db data/sit/sit_database.tsv || true
spoligotyper -i data/reads -o results/reads -t "$threads" -j 2 --operator validation --sit-db data/sit/sit_database.tsv || true
spoligotyper -i data/la_reads -o results/la_reads -t "$threads" -j 2 --operator validation --sit-db data/sit/sit_database.tsv || true
spoligotyper -i data/l1_reads -o results/l1_reads -t "$threads" -j 2 --operator validation --sit-db data/sit/sit_database.tsv || true
python3 check_results.py > results/validation.md
cat results/validation.md
