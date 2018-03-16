#!/usr/bin/python
import pysam
import re
import sys
import time
import os

from classes import read as r
from classes import segment as s

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

coverages = []
reads = {}
segments = {}
segmentID = 1

def parse_bam():
    """
    Reads bam file and saves reads and their segments in objects of the Read en Segment classes.
    :param bamfile used to open bam file:
    """
    global sample_name, header, segmentID, bam
    sys.stderr.write(time.strftime("%c") + " Busy with parsing bam file...\n")
    bam = pysam.AlignmentFile(NanoSV.opts_bam, 'rb')
    header = bam.header
    print( header )
    if 'HD' in header:
        if not header['HD']['SO'] == 'coordinate':
            sys.exit('The bam file is not coordinate sorted')
    if 'RG' in header:
        if type(header['RG']) is list:
            sample_name = header['RG'][0]['SM']
        else:
            sample_name = header['RG']['SM']
    else:
        sample_name = re.sub('(\.sorted)?\.bam$', '', str(NanoSV.opts_bam))

    for line in bam:
        if line.query_name in reads:
            read = reads[line.query_name]
        else:
            read = r.Read(line.query_name, line.infer_read_length())
            reads[line.query_name] = read
        
        if line.flag & 4 or line.mapping_quality < NanoSV.opts_min_mapq:
            continue
        segment = s.Segment(segmentID, line.query_name, line.flag, line.reference_name, line.reference_start+1, line.mapping_quality,
                            line.query_alignment_length)
        segment.end = line.reference_start + line.reference_length
        if line.has_tag('MD'):
            matches = sum( map(int, re.findall(r"(\d+)", line.get_tag('MD') )) )
            segment.pid = format(matches / segment.length, '.3f')
        else:
            segment.pid = format(line.get_cigar_stats()[0][7] / segment.length, '.3f')
            if segment.pid == "0.000":
                segment.pid = format(line.get_cigar_stats()[0][0] / segment.length, '.3f')
        if line.flag & 16:
            if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
                segment.clip = line.cigartuples[-1][1]
            else:
                segment.clip = 0
            if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
                segment.clip_2 = line.cigartuples[0][1]
            else:
                segment.clip_2 = 0
        else:
            if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
                segment.clip = line.cigartuples[0][1]
            else:
                segment.clip = 0
            if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
                segment.clip_2 = line.cigartuples[-1][1]
            else:
                segment.clip_2 = 0
        if float(segment.pid) < NanoSV.opts_min_pid:
            continue
        read.addSegment(segment)
        segments[segmentID] = segment
        segmentID += 1

        

