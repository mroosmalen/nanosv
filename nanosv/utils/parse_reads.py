#!/usr/bin/python
import sys
import time
import os

from collections import defaultdict

from classes import breakpoint as b
from utils import parse_bam as bam

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

breakpoints = {}
breakpoints_region = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
hanging_breakpoints_region = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))

def parse_reads():
    """
    Loops through all reads and saves breakpoints in objects of the Breakpoint class
    """
    sys.stderr.write(time.strftime("%c") + " Busy with parsing read segments...\n")
    breakpointID = 1
    hanging_breakpointID = -1
    for qname in bam.reads:
        clips = sorted(bam.reads[qname].segments)
        if len(clips) == 1 or len(clips) > NanoSV.opts_max_split:
            continue
        for i in range(0, (len(clips) - 1)):
            i2 = i + 1

            segment_1 = bam.segments[bam.reads[qname].segments[clips[i]]]
            segment_2 = bam.segments[bam.reads[qname].segments[clips[i2]]]
            segment_1.setPlength(bam.reads[qname].length)
            segment_2.setPlength(bam.reads[qname].length)
            breakpoint = b.Breakpoint(breakpointID, segment_1, segment_2)
            breakpointID += 1
            
            gap = (clips[i2] - (clips[i] + segment_1.length))
            breakpoint.setGap(gap)
            breakpoint.setBreakpoint(segment_1, segment_2)
            if segment_1.rname > segment_2.rname:
                breakpoint.switchSegments()
            elif segment_1.rname == segment_2.rname and breakpoint.segment_1["pos"] > breakpoint.segment_2["pos"]:
                breakpoint.switchSegments()

            breakpoint.setSVtype()
            breakpoints[breakpoint.id] = breakpoint
            values = ( breakpoint.svtype, str(breakpoint.segment_1["rname"]), str(breakpoint.segment_2["rname"]), str(breakpoint.segment_1["flag"]), str(breakpoint.segment_2["flag"]) )
            breakpoints_region["\t".join(values)][breakpoint.segment_1["pos"]][breakpoint.id] = 1
            
            if i == 0 and segment_1.clip >= NanoSV.opts_hanging_length:
                hanging_breakpoint_pos = segment_1.pos
                if segment_1.flag & 16:
                    hanging_breakpoint_pos = segment_1.end
                    hanging_breakpoints_region[segment_1.rname]['T'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_1.id
                    hanging_breakpointID -= 1
                else:
                    hanging_breakpoints_region[segment_1.rname]['H'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_1.id
                    hanging_breakpointID = hanging_breakpointID - 1
            if i2 == len(clips) - 1 and segment_2.clip_2 >= NanoSV.opts_hanging_length:
                hanging_breakpoint_pos = segment_2.end
                if segment_2.flag & 16:
                    hanging_breakpoint_pos = segment_2.pos
                    hanging_breakpoints_region[segment_2.rname]['H'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_2.id
                    hanging_breakpointID -= 1
                else:
                    hanging_breakpoints_region[segment_2.rname]['T'][hanging_breakpoint_pos][
                        hanging_breakpointID] = segment_2.id
                    hanging_breakpointID -= 1
