#!/usr/bin/python
import sys
import time
import os

from utils import parse_bam as bam
from utils import parse_breakpoints as breakpoint
from utils import phasing as phasing
from utils import create_vcf as c_vcf

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

def parse_svs():
    """
    Finds reference reads and prepares to print SVs to VCF
    """

    sys.stderr.write(time.strftime("%c") + " Busy with reference reads...\n")
    prev_rname = 0
    for chromosome in sorted(bam.segments):
        for position in sorted(bam.segments[chromosome]):
            for qname_clip in sorted(bam.segments[chromosome][position]):
                if bam.segments[chromosome][position][qname_clip].rname != prev_rname and prev_rname != 0:
                    del breakpoint.structural_variants_region[prev_rname]
                for sv_min in sorted(breakpoint.structural_variants_region[bam.segments[chromosome][position][qname_clip].rname]):
                    if sv_min <= bam.segments[chromosome][position][qname_clip].pos:
                        del breakpoint.structural_variants_region[bam.segments[chromosome][position][qname_clip].rname][sv_min]
                        continue
                    if sv_min > bam.segments[chromosome][position][qname_clip].end:
                        break
                    for sv_max in sorted(breakpoint.structural_variants_region[bam.segments[chromosome][position][qname_clip].rname][sv_min]):
                        if sv_max >= bam.segments[chromosome][position][qname_clip].end:
                            break
                        for sv_id in breakpoint.structural_variants_region[bam.segments[chromosome][position][qname_clip].rname][sv_min][sv_max]:
                            x = breakpoint.structural_variants_region[bam.segments[chromosome][position][qname_clip].rname][sv_min][sv_max][sv_id]
                            if sv_min > (bam.segments[chromosome][position][qname_clip].pos + NanoSV.opts_refread_distance) and \
                                    sv_max < (bam.segments[chromosome][position][qname_clip].end - NanoSV.opts_refread_distance):
                                breakpoint.structural_variants[sv_id].ref_qname_clips[x].append([chromosome, position, qname_clip])
                                breakpoint.structural_variants[sv_id].format['DR'][x] += 1
                                breakpoint.structural_variants[sv_id].format['RO'][x] += (1 - 10 ** (-bam.segments[chromosome][position][qname_clip].mapq / 10.0))
                prev_rname = bam.segments[chromosome][position][qname_clip].rname

    sys.stderr.write(time.strftime("%c") + " Busy with printing to vcf...\n")

    for sv_id in sorted(breakpoint.structural_variants):
        if breakpoint.structural_variants[sv_id].set != 1:
            breakpoint.structural_variants[sv_id].setArguments(NanoSV.opts_depth_support)
        if breakpoint.structural_variants[sv_id].format['GT'] == '0/0':
            continue
        for sv_id_2 in sorted(breakpoint.structural_variants_region_2[breakpoint.structural_variants[sv_id].chr]):
            if sv_id_2 <= sv_id:
                continue
            if breakpoint.structural_variants[sv_id_2].set != 1:
                breakpoint.structural_variants[sv_id_2].setArguments(NanoSV.opts_depth_support)
            if breakpoint.structural_variants[sv_id_2].format['GT'] == '0/0':
                continue
            if breakpoint.structural_variants[sv_id].chr != breakpoint.structural_variants[sv_id_2].chr:
                continue
            if breakpoint.structural_variants[sv_id].chr2 != breakpoint.structural_variants[sv_id_2].chr2:
                continue
            if abs(breakpoint.structural_variants[sv_id].pos - breakpoint.structural_variants[sv_id_2].pos) <= NanoSV.opts_window_size and abs(breakpoint.structural_variants[sv_id].info['END'] - breakpoint.structural_variants[sv_id_2].info['END']) <= NanoSV.opts_window_size:
                breakpoint.structural_variants[sv_id].SVcluster += 1
                breakpoint.structural_variants[sv_id_2].SVcluster += 1
            if breakpoint.structural_variants[sv_id_2].info['SVTYPE'] != "BND":
                continue
            if breakpoint.structural_variants[sv_id].info['SVTYPE'] != "BND":
                continue
            if (not breakpoint.structural_variants[sv_id].flag1 & 16 and not breakpoint.structural_variants[sv_id].flag2 & 16) and not (breakpoint.structural_variants[sv_id_2].flag1 & 16 and breakpoint.structural_variants[sv_id_2].flag2 & 16):
                continue
            if (not breakpoint.structural_variants[sv_id].flag1 & 16 and breakpoint.structural_variants[sv_id].flag2 & 16) and not (breakpoint.structural_variants[sv_id_2].flag1 & 16 and not breakpoint.structural_variants[sv_id_2].flag2 & 16):
                continue
            if (breakpoint.structural_variants[sv_id].flag1 & 16 and not breakpoint.structural_variants[sv_id].flag2 & 16) and not (not breakpoint.structural_variants[sv_id_2].flag1 & 16 and breakpoint.structural_variants[sv_id_2].flag2 & 16):
                continue
            if (breakpoint.structural_variants[sv_id].flag1 & 16 and breakpoint.structural_variants[sv_id].flag2 & 16) and not (not breakpoint.structural_variants[sv_id_2].flag1 & 16 and not breakpoint.structural_variants[sv_id_2].flag2 & 16):
                continue
            if abs(breakpoint.structural_variants[sv_id].pos - breakpoint.structural_variants[sv_id_2].pos) > NanoSV.opts_mate_distance:
                continue
            if abs(breakpoint.structural_variants[sv_id].info['END'] - breakpoint.structural_variants[sv_id_2].info['END']) > NanoSV.opts_mate_distance:
                continue
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

        windows = []
        if breakpoint.structural_variants[sv_id].flag1 == 16:
            windows.append([breakpoint.structural_variants[sv_id].pos, breakpoint.structural_variants[sv_id].pos + NanoSV.opts_phasing_window])
        elif breakpoint.structural_variants[sv_id].flag1 == 0:
            windows.append([breakpoint.structural_variants[sv_id].pos - NanoSV.opts_phasing_window, breakpoint.structural_variants[sv_id].pos])
        if breakpoint.structural_variants[sv_id].flag2 == 0:
            windows.append([breakpoint.structural_variants[sv_id].info['END'], breakpoint.structural_variants[sv_id].info['END'] + NanoSV.opts_phasing_window])
        elif breakpoint.structural_variants[sv_id].flag2 == 16:
            windows.append([breakpoint.structural_variants[sv_id].info['END'] - NanoSV.opts_phasing_window, breakpoint.structural_variants[sv_id].info['END']])
        if breakpoint.structural_variants[sv_id].format['GT'] == '0/1' and NanoSV.opts_phasing_on:
            phasing_result = phasing.make_matrix(sv_id, windows)
            breakpoint.structural_variants[sv_id].info['PURITY_SCORE'] = phasing_result[0]
            breakpoint.structural_variants[sv_id].info['PHASING_SCORE'] = phasing_result[1]
            breakpoint.structural_variants[sv_id].info['SNPS_USED'] = phasing_result[3]
            breakpoint.structural_variants[sv_id].info['PHASING_PVALUE'] = phasing_result[2]
        breakpoint.structural_variants[sv_id].printVCF(c_vcf.vcf_writer)
