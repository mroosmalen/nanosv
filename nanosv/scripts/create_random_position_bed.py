#!/usr/bin/python
# Create random position bed for hg19

import random

pick_random = 100000
mask_regions = dict()

genome = {
    1: 249250621,
    2: 243199373,
    3: 198022430,
    4: 191154276,
    5: 180915260,
    6: 171115067,
    7: 159138663,
    8: 146364022,
    9: 141213431,
    10: 135534747,
    11: 135006516,
    12: 133851895,
    13: 115169878,
    14: 107349540,
    15: 102531392,
    16: 90354753,
    17: 81195210,
    18: 78077248,
    19: 59128983,
    20: 63025520,
    21: 48129895,
    22: 51304566,
    23: 155270560
}

genome_length = sum(genome.values())

simple_repeats_file = '/hpc/cog_bioinf/kloosterman/users/mroosmalen/data/simpleRepeats.bed'
gaps_file = '/hpc/cog_bioinf/kloosterman/users/mroosmalen/data/gap.bed'


def read_bed(file):
    with (open(file, 'r')) as bed:
        for line in bed:
            line = line.rstrip()
            chr, start, end = line.split("\t")
            chr = chr.replace('chr', '')
            if chr == 'X':
                chr = 23
            try:
                chr = int(chr)
            except:
                continue
            if chr not in genome:
                continue
            if chr not in mask_regions:
                mask_regions[chr] = dict()
            if start not in mask_regions[chr]:
                mask_regions[chr][int(start)] = dict()
            mask_regions[chr][int(start)][int(end)] = 1


read_bed(simple_repeats_file)
read_bed(gaps_file)

random_positions = dict()

while len(random_positions.keys()) < pick_random:
    randnum = random.randrange(1, genome_length)
    tmp_length = 0
    for chr, chrlength in sorted(genome.iteritems()):
        if randnum < tmp_length + chrlength:
            randchr = chr
            randpos = randnum - tmp_length
            break
        tmp_length += chrlength

    mask = False
    for mask_start in sorted(mask_regions[randchr]):
        if randpos < mask_start:
            break
        for mask_end in sorted(mask_regions[randchr][mask_start]):
            if randpos <= mask_end:
                mask = True
                break
    if not mask:
        random_positions["\t".join(
            [str(randchr), str(randpos),
             str(randpos + 1)])] = 1

print("\n".join(random_positions.keys()))
