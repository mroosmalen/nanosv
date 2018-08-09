#!/usr/bin/python
# Create random position bed for hg19

import random

pick_random = 100000
mask_regions = {}

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
            ch, start, end = line.split("\t")
            ch = ch.replace('chr', '')
            if ch == 'X':
                ch = 23
            elif ch == 'Y':
                ch = 24
            try:
                ch = int(ch)
            except:
                continue
            if ch not in genome:
                continue
            if ch not in mask_regions:
                mask_regions[ch] = {}
            if start not in mask_regions[ch]:
                mask_regions[ch][int(start)] = {}
            mask_regions[ch][int(start)][int(end)] = 1


read_bed(simple_repeats_file)
read_bed(gaps_file)

random_positions = {}

while len(random_positions.keys()) < pick_random:
    randnum = random.randrange(1, genome_length)
    tmp_length = 0
    for ch, chrlength in sorted(genome.iteritems()):
        if randnum < tmp_length + chrlength:
            randchr = ch
            randchr = randchr.replace('chr','')
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
        if randchr == '23':
            ranchr = 'X'
        elif ranchr == '24':
            randchr = 'Y'
        random_positions["\t".join(
            [str(randchr), str(randpos),
             str(randpos + 1)])] = 1

print("\n".join(random_positions.keys()))
