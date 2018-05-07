#!/usr/bin/env bash
# Create random position bed for mm10

set -euo pipefail

die () {
    echo >&2 "$@"
    exit 1
}

[[ "$#" -eq 1 ]] || die "1 argument required (path to .genome), $# provided"
[[ -f "$1" ]] || die "File $1 does not exist"

if [[ ! -f gap.bed ]]
then
    mysql \
        --user=genome \
        --host=genome-mysql.soe.ucsc.edu \
        -A \
        -P 3306 \
        -B \
        -N \
        -e "SELECT chrom, chromStart, chromEnd from simpleRepeat;" \
        mm10 > simpleRepeats.bed
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
        mm10 > gap.bed
fi

# Combine gaps and simple repeats to single exclusion bed
cat gap.bed simpleRepeats.bed > exclude.bed

# Choose 1M random spots in the genome
bedtools random \
    -l 1 \
    -g "$1" > raw_genome_sample.bed

# Shuffle positions while excluding gaps/simple repeats
# Should not have duplicates
bedtools shuffle \
    -excl exclude.bed \
    -noOverlapping \
    -i raw_genome_sample.bed \
    -g "$1" > mm10_genome_sample.bed
