#!/usr/bin/env bash
# Create random position bed for any genome
# First argument should UCSC genome name (e.g. mm10 or hg38)
# Second argument is Bedtools genome file [chr  length]



set -euo pipefail

die () {
    echo >&2 "$@"
    exit 1
}

[[ "$#" -eq 2 ]] || die "2 arguments required (Genome name and path to .genome), $# provided"
[[ -f "$2" ]] || die "File $2 does not exist"

if [[ ! -f simpleRepeats.bed ]]
then
    mysql \
        --user=genome \
        --host=genome-mysql.soe.ucsc.edu \
        -A \
        -P 3306 \
        -B \
        -N \
        -e "SELECT chrom, chromStart, chromEnd from simpleRepeat;" \
        "$1" > simpleRepeats.bed
fi

if [[ ! -f gap.bed ]]
then
    mysql \
        --user=genome \
        --host=genome-mysql.soe.ucsc.edu \
        -A \
        -P 3306 \
        -B \
        -N \
        -e "SELECT chrom, chromStart, chromEnd from gap;" \
        "$1" > gap.bed
fi

# Combine gaps and simple repeats to single exclusion bed
cat gap.bed simpleRepeats.bed > exclude.bed

# Choose 1M random spots in the genome
bedtools random \
    -l 1 \
    -g "$2" > raw_genome_sample.bed

# Shuffle positions while excluding gaps/simple repeats
# Should not have duplicates
bedtools shuffle \
    -excl exclude.bed \
    -noOverlapping \
    -i raw_genome_sample.bed \
    -g "$2" > "$1"_genome_sample.bed
