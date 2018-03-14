#!/usr/bin/python
import re
import sys
import time
import os

from collections import defaultdict

from classes import sv as svclass
from utils import parse_reads as read
from utils import parse_bam as bam

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

svID = 1
structural_variants = {}
structural_variants_region = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
structural_variants_region_2 = defaultdict(lambda: defaultdict(int))


def addSVInfo(sv):
    """
    Adds info to the newly made SV object
    :param sv is SV-object:
    """
    if sum(sv.format['DV']) >= NanoSV.opts_cluster_count * 2:
        flag_list = [sv.flag1, sv.flag2]
        for flag_idx in range(len(flag_list)):            
            if ( flag_idx == 0 and not flag_list[flag_idx] & 16 ) or ( flag_idx == 1 and flag_list[flag_idx] & 16 ):
                head_tail = 'T'
            else:
                head_tail = 'H'

            for hanging_pos in sorted(read.hanging_breakpoints_region[sv.chr][head_tail]):
                if hanging_pos < (min(sv.pos) - NanoSV.opts_cluster_distance):
                    continue
                if hanging_pos > (max(sv.pos) + NanoSV.opts_cluster_distance):
                    break
                for hanging_id in read.hanging_breakpoints_region[sv.chr][head_tail][hanging_pos]:
                    segment_id = read.hanging_breakpoints_region[sv.chr][head_tail][hanging_pos][hanging_id]
                    if segment_id == 0: continue
                    sv.format['HR'][flag_idx] += 1
                    sv.format['VO'][flag_idx] += (1 - 10 ** (-bam.segments[segment_id].mapq / 10.0))
                    sv.addInfoField("PID", [[bam.segments[segment_id].pid], None])
                    sv.addInfoField("MAPQ", [[bam.segments[segment_id].mapq], None])
                    sv.addInfoField("PLENGTH", [[bam.segments[segment_id].plength], None])
                    sv.addInfoField("RLENGTH", [bam.reads[bam.segments[segment_id].qname].length])
                    if re.match("/2D_2d$/", bam.segments[segment_id].qname):
                        sv.addInfoField("RT", ["2d"])
                    elif re.match("/2D_complement$/", bam.segments[segment_id].qname):
                        sv.addInfoField("RT", ["complement"])
                    else:
                        sv.addInfoField("RT", ["template"])
        structural_variants[sv.id] = sv
        structural_variants_region[sv.chr][min(sv.pos)][max(sv.pos)][sv.id] = 0
        structural_variants_region[sv.chr2][min(sv.info['END'])][max(sv.info['END'])][sv.id] = 1
        structural_variants_region_2[sv.chr][sv.id] = 1


def parse_breakpoints_2(breakpoints_region_2):
    """
    Creates an SV object for every SV.
    :param breakpoints_region_2:
    """
    global svID
    prev_pos_2 = -1
    for pos_2 in sorted(breakpoints_region_2):
        for pos_1 in sorted(breakpoints_region_2[pos_2]):
            for breakpoint_id in breakpoints_region_2[pos_2][pos_1]:
                breakpoint = read.breakpoints[breakpoint_id]
                if (prev_pos_2 == -1):
                    sv = svclass.SV(svID, breakpoint)
                    svID += 1
                elif abs(breakpoint.segment_2["pos"] - prev_pos_2) <= NanoSV.opts_cluster_distance:
                    sv.addBreakpoint(breakpoint)
                else:
                    addSVInfo(sv)
                    sv = svclass.SV(svID, breakpoint)
                    svID += 1

                sv.addInfoField("PID",[ [bam.segments[breakpoint.segment_1["id"]].pid],[bam.segments[breakpoint.segment_2["id"]].pid] ])
                sv.addInfoField("MAPQ",[ [bam.segments[breakpoint.segment_1["id"]].mapq],[bam.segments[breakpoint.segment_2["id"]].mapq] ])
                sv.addInfoField("PLENGTH",[ [bam.segments[breakpoint.segment_1["id"]].plength],[bam.segments[breakpoint.segment_2["id"]].plength] ])
                sv.addInfoField("RLENGTH",[ bam.reads[bam.segments[breakpoint.segment_2["id"]].qname].length ])
                sv.addInfoField("GAP", [breakpoint.gap])

                if re.match("/2D_2d$/", bam.segments[breakpoint.segment_2["id"]].qname):
                    sv.addInfoField("RT", ["2d"])
                elif re.match("/2D_complement$/", bam.segments[breakpoint.segment_2["id"]].qname):
                    sv.addInfoField("RT", ["complement"])
                else:
                    sv.addInfoField("RT", ["template"])
                prev_pos_2 = int(pos_2)
    addSVInfo(sv)


def parse_breakpoints():
    """
    Loops through breakpoints_region dict and calls parse_breakpoints_2
    """
    sys.stderr.write(time.strftime("%c") + " Busy with parsing breakpoints...\n")

    for region in read.breakpoints_region:
        prev_pos_1 = -1
        breakpoints_region_2 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        for pos_1 in sorted(read.breakpoints_region[region]):
            for breakpoint_id in read.breakpoints_region[region][pos_1]:
                if prev_pos_1 == -1:
                    breakpoints_region_2[read.breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                elif abs(pos_1 - prev_pos_1) <= NanoSV.opts_cluster_distance:
                    breakpoints_region_2[read.breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                else:
                    parse_breakpoints_2(breakpoints_region_2)
                    breakpoints_region_2.clear()
                    breakpoints_region_2[read.breakpoints[breakpoint_id].segment_2["pos"]][pos_1][breakpoint_id] = 1
                prev_pos_1 = int(pos_1)
        parse_breakpoints_2(breakpoints_region_2)

    return structural_variants, structural_variants_region, structural_variants_region_2
