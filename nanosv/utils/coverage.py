#!/usr/bin/python

import os
import time
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

coverages = []

def calculate_coverage_bed():
    """
    Create list of coverages based on locations in bed file.
    """
    sys.stderr.write(time.strftime("%c") + " Busy with calculating the coverage distribution...\n")
    if 'sambamba' in NanoSV.opts_sambamba:
        cmd = NanoSV.opts_sambamba + " depth base --min-coverage=0 " + NanoSV.opts_bam + " -L " + NanoSV.opts_bed + " --nthreads=" + str(NanoSV.opts_threads) + " 2> /dev/null | awk '{if (NR!=1) print $3}'"
    elif 'samtools' in NanoSV.opts_sambamba:
        cmd = NanoSV.opts_sambamba + " depth " + NanoSV.opts_bam + " -b " + NanoSV.opts_bed + " | awk '{print $3}'"
    else:
        sys.exit('No sambamba or samtools found')

    with os.popen(cmd) as coverageOutput:
        for coverage in coverageOutput:
            if coverage != "" and coverage != "\n":
                coverages.append(int(coverage))

    if len(coverages) == 0:
        sys.exit("Can't calculate coverage distribution. The bed file may be inappropriate for your bam file.")
