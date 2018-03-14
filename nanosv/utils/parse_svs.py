#!/usr/bin/python
import re
import sys
import time
import os 
import pysam

from utils import parse_bam as bam
from utils import parse_reads as read
from utils import parse_breakpoints as breakpoint
from utils import create_vcf as vcf

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

def parse_svs():
    """
    Finds reference reads and prepares to print SVs to VCF
    """
    sys.stderr.write(time.strftime("%c") + " Busy with reference reads...\n")
    prev_rname = 0
    for segment_id in sorted(bam.segments):
        if bam.segments[segment_id].rname != prev_rname and prev_rname != 0:
            del breakpoint.structural_variants_region[prev_rname]
        for sv_min in sorted(breakpoint.structural_variants_region[bam.segments[segment_id].rname]):
            if sv_min <= bam.segments[segment_id].pos:
                del breakpoint.structural_variants_region[bam.segments[segment_id].rname][sv_min]
                continue
            if sv_min > bam.segments[segment_id].end:
                break
            for sv_max in sorted(breakpoint.structural_variants_region[bam.segments[segment_id].rname][sv_min]):
                if sv_max >= bam.segments[segment_id].end:
                    break
                for sv_id in breakpoint.structural_variants_region[bam.segments[segment_id].rname][sv_min][sv_max]:
                    x = breakpoint.structural_variants_region[bam.segments[segment_id].rname][sv_min][sv_max][sv_id]
                    if sv_min > (bam.segments[segment_id].pos + NanoSV.opts_refread_distance) and \
                            sv_max < (bam.segments[segment_id].end - NanoSV.opts_refread_distance):
                        breakpoint.structural_variants[sv_id].format['DR'][x] += 1
                        breakpoint.structural_variants[sv_id].format['RO'][x] += (1 - 10 ** (-bam.segments[segment_id].mapq / 10.0))
        prev_rname = bam.segments[segment_id].rname

    sys.stderr.write(time.strftime("%c") + " Busy with printing to vcf...\n")

    for sv_id in sorted(breakpoint.structural_variants):
        if breakpoint.structural_variants[sv_id].set != 1: breakpoint.structural_variants[sv_id].setArguments(NanoSV.opts_depth_support)
        if breakpoint.structural_variants[sv_id].format['GT'] == '0/0': continue
        for sv_id_2 in sorted(breakpoint.structural_variants_region_2[breakpoint.structural_variants[sv_id].chr]):
            if sv_id_2 <= sv_id: continue
            if breakpoint.structural_variants[sv_id_2].set != 1: breakpoint.structural_variants[sv_id_2].setArguments(NanoSV.opts_depth_support)
            if breakpoint.structural_variants[sv_id_2].format['GT'] == '0/0': continue
            if breakpoint.structural_variants[sv_id].chr != breakpoint.structural_variants[sv_id_2].chr: continue
            if breakpoint.structural_variants[sv_id].chr2 != breakpoint.structural_variants[sv_id_2].chr2: continue
            if abs(breakpoint.structural_variants[sv_id].pos - breakpoint.structural_variants[sv_id_2].pos) <= NanoSV.opts_window_size and abs(breakpoint.structural_variants[sv_id].info['END'] - breakpoint.structural_variants[sv_id_2].info['END']) <= NanoSV.opts_window_size:
                breakpoint.structural_variants[sv_id].SVcluster += 1
                breakpoint.structural_variants[sv_id_2].SVcluster += 1
            if breakpoint.structural_variants[sv_id_2].info['SVTYPE'] != "BND": continue
            if breakpoint.structural_variants[sv_id].info['SVTYPE'] != "BND": continue
            if ( not breakpoint.structural_variants[sv_id].flag1 & 16 and not breakpoint.structural_variants[sv_id].flag2 & 16) and not (breakpoint.structural_variants[sv_id_2].flag1 & 16 and breakpoint.structural_variants[sv_id_2].flag2 & 16): continue
            if ( not breakpoint.structural_variants[sv_id].flag1 & 16 and breakpoint.structural_variants[sv_id].flag2 & 16) and not (breakpoint.structural_variants[sv_id_2].flag1 & 16 and not breakpoint.structural_variants[sv_id_2].flag2 & 16): continue
            if (breakpoint.structural_variants[sv_id].flag1 & 16 and not breakpoint.structural_variants[sv_id].flag2 & 16) and not ( not breakpoint.structural_variants[sv_id_2].flag1 & 16 and breakpoint.structural_variants[sv_id_2].flag2 & 16): continue
            if (breakpoint.structural_variants[sv_id].flag1 & 16 and breakpoint.structural_variants[sv_id].flag2 & 16) and not ( not breakpoint.structural_variants[sv_id_2].flag1 & 16 and not breakpoint.structural_variants[sv_id_2].flag2 & 16): continue
            if abs(breakpoint.structural_variants[sv_id].pos - breakpoint.structural_variants[sv_id_2].pos) > NanoSV.opts_mate_distance: continue
            if abs(breakpoint.structural_variants[sv_id].info['END'] - breakpoint.structural_variants[sv_id_2].info['END']) > NanoSV.opts_mate_distance: continue
            breakpoint.structural_variants[sv_id].info['MATEID'] = sv_id_2
            breakpoint.structural_variants[sv_id_2].info['MATEID'] = sv_id
        if breakpoint.structural_variants[sv_id].qual < NanoSV.opts_qual_flag:
            breakpoint.structural_variants[sv_id].filter.append("LowQual")
        if breakpoint.structural_variants[sv_id].SVcluster > NanoSV.opts_svcluster:
            breakpoint.structural_variants[sv_id].filter.append("SVcluster")
        if breakpoint.structural_variants[sv_id].info['GAP'] > NanoSV.opts_gap_flag and breakpoint.structural_variants[sv_id].info['SVTYPE'] != "INS":
            breakpoint.structural_variants[sv_id].filter.append("GAP")
        if breakpoint.structural_variants[sv_id].info['MAPQ'][0] < NanoSV.opts_mapq_flag or breakpoint.structural_variants[sv_id].info['MAPQ'][1] < NanoSV.opts_mapq_flag:
            breakpoint.structural_variants[sv_id].filter.append("MapQual")
        if breakpoint.structural_variants[sv_id].info['PID'][0] < NanoSV.opts_pid_flag or breakpoint.structural_variants[sv_id].info['PID'][1] < NanoSV.opts_pid_flag:
            breakpoint.structural_variants[sv_id].filter.append("PID")
        if abs(breakpoint.structural_variants[sv_id].info['CIPOS'][0])+breakpoint.structural_variants[sv_id].info['CIPOS'][1] > NanoSV.opts_ci_flag:
            breakpoint.structural_variants[sv_id].filter.append("CIPOS")
        if abs(breakpoint.structural_variants[sv_id].info['CIEND'][0])+breakpoint.structural_variants[sv_id].info['CIEND'][1] > NanoSV.opts_ci_flag:
            breakpoint.structural_variants[sv_id].filter.append("CIEND")
        if breakpoint.structural_variants[sv_id].info['SVTYPE'] == 'INS':
            breakpoint.structural_variants[sv_id].info['SVLEN'] = breakpoint.structural_variants[sv_id].info['GAP']
        breakpoint.structural_variants[sv_id].printVCF(vcf.vcf_writer)
