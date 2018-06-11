#!/usr/bin/python
import pysam
import re
import sys
import time
import os
import NanoSV

from classes import read as r
from classes import segment as s
from classes import variant as v

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))


coverages = []
reads = {}
segments = {}
tmp_variants = {}
variants = {}
segments_to_check = {}


def parse_bam():
    """
    Reads bam file and saves reads and their segments in objects of the Read en Segment classes.
    :param bamfile used to open bam file:
    """
    global sample_name, header, segmentID, bam, tmp_variants, segments_to_check
    segmentID = 1
    sys.stderr.write(time.strftime("%c") + " Busy with parsing bam file...\n")
    bam = pysam.AlignmentFile(NanoSV.opts_bam, 'rb')
    header = bam.header
    if not bam.has_index():
        sys.exit('The bam has no index file')
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
    for item in bam.header["SQ"]:
        variants[item['SN']] = {}
        for bin in range(0,int(item['LN']/NanoSV.opts_variant_bin_size)):
            variants[item['SN']][bin] = {}
    previous_refname = -1
    previous_cursor = -1
    for line in bam:
        if line.flag & 4:
            continue
        if segmentID % 100 == 0:
            sys.stderr.write(time.strftime("%c") + " " + str(segmentID) + " reads geladen\n")
        remove = [qname_clip for qname_clip in segments_to_check if line.reference_start > segments_to_check[qname_clip][1]]
        for q_clip in remove: del segments_to_check[q_clip]
        if line.query_name in reads:
            read = reads[line.query_name]
        else:
            read = r.Read(line.query_name, line.infer_read_length())
            reads[line.query_name] = read
        clip, clip_2 = calculate_clip(line)
        created_subsegments = search_for_indels(line, clip, clip_2)
        for segment in created_subsegments:
            if segment.rname != previous_refname:
                if NanoSV.opts_phasing_on:
                    tmp_variants = {}
                segments[segment.rname] = {}
            if not int(segment.pos) in segments[segment.rname]:
                segments[segment.rname][segment.pos] = {}
            read.addSegment(segment)
            segments[segment.rname][segment.pos][segment.id[-1]] = segment
            segments_to_check[segment.id[-1]] = [segment.pos, segment.end]
            previous_refname = line.reference_name
        if NanoSV.opts_phasing_on and line.seq and line.query_qualities:
            if line.has_tag('MD'):
                if not NanoSV.opts_snp_file:
                    find_variations_md(line.cigartuples, line.reference_name, line.reference_start + 1, line.query_qualities, line.seq, line.get_tag('MD'), created_subsegments)
            else:
                if not NanoSV.opts_snp_file:
                    find_variations_cigar(line.cigartuples, line.reference_name, line.reference_start + 1, line.query_qualities, line.seq, created_subsegments)
        if previous_cursor == -1:
            previous_cursor = line.reference_start
        if NanoSV.opts_phasing_on and line.seq and line.query_qualities and line.reference_start > previous_cursor:
            remove_variations(previous_cursor, line.reference_start, line.reference_name)
        previous_cursor = line.reference_start
    if NanoSV.opts_phasing_on and NanoSV.opts_snp_file:
        read_snp_vcf()
    write_bed()


def calculate_pid(line, query_alignment_length):
    """
    calculates percentage identity of alignment of read
    :param line: 
    :param query_alignment_length: 
    :return pid: 
    """
    if line.has_tag('MD'):
        matches = sum(map(int, re.findall(r"(\d+)", line.get_tag('MD'))))
        pid = format(matches / query_alignment_length, '.3f')
    else:
        pid = format(line.get_cigar_stats()[0][7] / query_alignment_length, '.3f')
        if pid == "0.000":
            pid = format(line.get_cigar_stats()[0][0] / query_alignment_length, '.3f')
    return pid


def calculate_clip(line):
    """
    calculates clip of segment
    :param line: 
    :return list with clip of both sides of the segment: 
    """
    if line.flag & 16:
        if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
            clip = line.cigartuples[-1][1]
        else:
            clip = 0
        if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
            clip_2 = line.cigartuples[0][1]
        else:
            clip_2 = 0
    else:
        if line.cigartuples[0][0] == 5 or line.cigartuples[0][0] == 4:
            clip = line.cigartuples[0][1]
        else:
            clip = 0
        if line.cigartuples[-1][0] == 5 or line.cigartuples[-1][0] == 4:
            clip_2 = line.cigartuples[-1][1]
        else:
            clip_2 = 0
    return [clip, clip_2]


def create_pattern():
    """
    creates a pattern to search for deletions and indels in cigar string using regular expressions
    :return pattern to search for with re.finditer(): 
    """
    pattern_list = [''] * len(NanoSV.opts_min_indel_size)
    for i in range(0, len(NanoSV.opts_min_indel_size)):
        j = len(NanoSV.opts_min_indel_size) - i - 1
        if j > 0:
            pattern_list[i] += '[' + str(int(NanoSV.opts_min_indel_size[i]) + 1) + '-9]\d{' + str(j) + ',}'
        for k in range(i + 1, len(pattern_list)):
            pattern_list[k] += '[' + str(int(NanoSV.opts_min_indel_size[i])) + '-9]'
    pattern_list[-1] += '[' + str(int(NanoSV.opts_min_indel_size[-1])) + '-9]'
    pattern_list.append("\d{" + str(len(NanoSV.opts_min_indel_size) + 1) + ",}")
    pattern = re.compile(r'' + "(" + "|".join(pattern_list) + ")([DI])" + '')
    return pattern


def search_for_indels(line, clip, clip_2):
    """
    Looks for indels that could be broken up to create breakpoints. Breakpoints are detected by NanoSV
    (Gapped alignment fix)
    :param line: 
    :param clip: 
    :param clip_2: 
    :return list of created subsegments: 
    """
    cigarstring = line.cigarstring
    reference_start = line.reference_start + 1
    line_query_alignment_length = line.query_alignment_length
    query_alignment_length = 0
    c_start = 0
    pattern = create_pattern()
    i = 0
    created_subsegments = []
    for m in re.finditer(pattern, line.cigarstring):
        line.cigarstring = cigarstring[c_start:m.start(0)]
        reference_end = reference_start + line.get_cigar_stats()[0][0] + line.get_cigar_stats()[0][2] + \
                        line.get_cigar_stats()[0][7] + line.get_cigar_stats()[0][8] - 1
        if line.flag & 16:
            clip_2 = clip_2 + query_alignment_length
            query_alignment_length = line.get_cigar_stats()[0][0] + line.get_cigar_stats()[0][1] + \
                                     line.get_cigar_stats()[0][7] + line.get_cigar_stats()[0][8]
            if i == 0:
                clip = clip + line_query_alignment_length - query_alignment_length
            else:
                clip = clip - query_alignment_length
            i += 1
        else:
            clip = clip + query_alignment_length

            query_alignment_length = line.get_cigar_stats()[0][0] + line.get_cigar_stats()[0][1] + \
                                     line.get_cigar_stats()[0][7] + line.get_cigar_stats()[0][8]
            if i == 0:
                clip_2 = clip_2 + line_query_alignment_length - query_alignment_length
            else:
                clip_2 = clip_2 - query_alignment_length
            i += 1
        if query_alignment_length == 0:
            reference_start = int(reference_end) + 1
            if line.flag & 16:
                if m.group(2) == "I":
                    clip -= int(m.group(1))
            else:
                if m.group(2) == "I":
                    clip += int(m.group(1))
            if m.group(2) == "D":
                reference_start += int(m.group(1))
            c_start = m.end(0)
            continue
        pid = keep_segment(line, query_alignment_length)
        if not pid:
            reference_start = int(reference_end) + 1
            if line.flag & 16:
                if m.group(2) == "I":
                    clip -= int(m.group(1))
            else:
                if m.group(2) == "I":
                    clip += int(m.group(1))
            if m.group(2) == "D":
                reference_start += int(m.group(1))
            c_start = m.end(0)
            continue
        segment = s.Segment(line.query_name, line.flag, line.reference_name, reference_start, line.mapping_quality, query_alignment_length, clip, clip_2, pid)
        segment.end = reference_end
        created_subsegments.append(segment)
        reference_start = int(reference_end) + 1
        if line.flag & 16:
            if m.group(2) == "I":
                clip -= int(m.group(1))
        else:
            if m.group(2) == "I":
                clip += int(m.group(1))
        if m.group(2) == "D":
            reference_start += int(m.group(1))
        c_start = m.end(0)
    line.cigarstring = cigarstring[c_start:]
    reference_end = reference_start + line.get_cigar_stats()[0][0] + line.get_cigar_stats()[0][2] + \
                    line.get_cigar_stats()[0][7] + line.get_cigar_stats()[0][8] - 1
    if line.flag & 16:
        if i != 0:
            clip_2 = clip_2 + query_alignment_length
            query_alignment_length = line.get_cigar_stats()[0][0] + line.get_cigar_stats()[0][1] + line.get_cigar_stats()[0][7] + line.get_cigar_stats()[0][8]
            clip = clip - query_alignment_length
        else:
            query_alignment_length = line_query_alignment_length
    else:
        if i != 0:
            clip = clip + query_alignment_length
            query_alignment_length = line.get_cigar_stats()[0][0] + line.get_cigar_stats()[0][1] + \
                                     line.get_cigar_stats()[0][7] + line.get_cigar_stats()[0][8]
            clip_2 = clip_2 - query_alignment_length
        else:
            query_alignment_length = line_query_alignment_length
    pid = keep_segment(line, query_alignment_length)
    if pid:
        segment = s.Segment(line.query_name, line.flag, line.reference_name, reference_start, line.mapping_quality, query_alignment_length, clip, clip_2, pid)
        segment.end = reference_end
        created_subsegments.append(segment)
    return created_subsegments


def read_snp_vcf():
    """Reads vcf/bedfile with pre-set SNP positions and calls find_SNPs() to save the variants"""
    with open(NanoSV.opts_snp_file, "r") as snp_file:
        snp_nummer = 1
        for line in snp_file:
            snp_nummer += 1
            line = line.rstrip()
            if not line.startswith("#"):
                columns = line.split("\t")
                if len(columns) > 1:
                    find_SNPs(columns[0], int(columns[1]))


def find_SNPs(chromosome, snp_position):
    """
    Looks for SNPs on given position in the genome using pileup.
    :param chromosome: 
    :param snp_position: 
    """
    base_ratios = {'A': [0, 0], 'C': [0, 0], 'G': [0, 0], 'T': [0, 0], '=': [0, 0]}
    deletions = 0
    total_n = 0
    variant = v.Variant(chromosome, int(snp_position))
    for pileupcolumn in bam.pileup(chromosome, int(snp_position)-1, int(snp_position), truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not keep_segment(pileupread.alignment):
                continue
            clip, clip_2 = calculate_clip(pileupread.alignment)
            if pileupread.is_del:
                variant.add_segment([pileupread.alignment.reference_name, pileupread.alignment.reference_start, str(pileupread.alignment.query_name) + ";" + str(clip)], '-')
                deletions += 1
                total_n += 1
            if not pileupread.is_del and not pileupread.is_refskip:
                if pileupread.alignment.query_qualities[pileupread.query_position] >= NanoSV.opts_min_base_qual_ph:
                    base_ratios[pileupread.alignment.query_sequence[pileupread.query_position]][0] += 1
                else:
                    base_ratios[pileupread.alignment.query_sequence[pileupread.query_position]][1] += 1
                variant.add_segment([pileupread.alignment.reference_name, pileupread.alignment.reference_start, str(pileupread.alignment.query_name) + ";" + str(clip)], [pileupread.alignment.query_sequence[pileupread.query_position], pileupread.alignment.query_qualities[pileupread.query_position]])
                total_n += 1
    if deletions < (NanoSV.opts_max_deletions * total_n):
        haplotypes = sorted(base_ratios.items(), key=lambda x: sum(x[1]))[-2:]
        try:
            if haplotypes[0][1][0] / sum(haplotypes[0][1]) > NanoSV.opts_min_occurences_of_highq_var and haplotypes[1][1][0] / sum(haplotypes[1][1]) > NanoSV.opts_min_occurences_of_highq_var:
                if sum(haplotypes[0][1]) / (sum(haplotypes[1][1]) + sum(haplotypes[0][1])) > NanoSV.opts_min_occurences_of_var:
                    bin = int(int(snp_position) / NanoSV.opts_variant_bin_size)
                    variants[chromosome][bin][int(snp_position)] = variant
        except ZeroDivisionError:
                ""


def keep_segment(line, query_alignment_length):
    """
    Calls calculate_pid() and checks if pid meets requirement to use read in NanoSV. If so, returns pid; else False
    :param line: 
    :param query_alignment_length: 
    :return pid: 
    """
    if line.mapping_quality < NanoSV.opts_min_mapq:
        return False
    pid = calculate_pid(line, query_alignment_length)
    if float(pid) < NanoSV.opts_min_pid:
        return False
    return pid


def find_variations_md(cigartuples, chromosome, pos, qual, seq, md, subsegments):
    """
    Loops through read to find variant base. This function is used with BWA-MEM, Minimap2 and NGMLR bam files
    :param cigartuples: 
    :param chromosome: 
    :param pos: 
    :param qual: 
    :param seq: 
    :param md: 
    :param subsegments: 
    """
    global tmp_variants
    md_tag = re.split("\^\D+", md)
    read_cursor = 0
    ref_cursor = 0
    md_cursor = 0
    md_string = ''
    for tuple in cigartuples:
        if tuple[0] == 4:
            read_cursor += tuple[1]
        if tuple[0] == 1:
            read_cursor += tuple[1]
        if tuple[0] == 2:
            md_cursor += 1
            md_string = ''
            for deletion in range(tuple[1]):
                if (pos + 1 + ref_cursor) not in tmp_variants:
                    tmp_variants[pos + 1 + ref_cursor] = v.Variant(chromosome, pos + 1 + ref_cursor)
                for segment in subsegments:
                    if segment.pos <= (pos + 1 + ref_cursor) <= segment.end:
                        tmp_variants[pos + 1 + ref_cursor].add_segment(segment.id, ['-'])
                        break
                    if segment.end < ref_cursor:
                        del subsegments[subsegments.index(segment)]
                ref_cursor += 1
        if tuple[0] == 0:
            if re.search("\D+", md_tag[md_cursor]):
                if md_string == '':
                    for m in re.split("\D", md_tag[md_cursor]):
                        if md_string != '':
                            md_string += 'X'
                        md_string += "=" * int(m)
                for seq_mismatch in re.finditer('X', md_string[:tuple[1]]):
                    if (pos + 1 + ref_cursor) not in tmp_variants:
                        tmp_variants[pos + 1 + ref_cursor] = v.Variant(chromosome, pos + 1 + ref_cursor)
                    for segment in subsegments:
                        if segment.pos <= (pos + 1 + ref_cursor) <= segment.end:
                            tmp_variants[pos + 1 + ref_cursor].add_segment(segment.id, [seq[read_cursor + seq_mismatch.start()], qual[read_cursor + seq_mismatch.start()]])
                            break
                        if segment.end < ref_cursor:
                            del subsegments[subsegments.index(segment)]
                md_string = md_string[tuple[1]:]
            ref_cursor += tuple[1]
            read_cursor += tuple[1]


def find_variations_cigar(cigartuples, chromosome, pos, qual, seq, subsegments):
    """
    Loops through read to find variant base. This function is used with LAST bam files
    :param cigartuples: 
    :param chromosome: 
    :param pos: 
    :param qual: 
    :param seq: 
    :param subsegments: 
    """
    global tmp_variants
    ref_cursor = (int(pos))
    read_cursor = 0
    for tuple in cigartuples:
        if tuple[0] == 4:
            read_cursor += tuple[1]
        if tuple[0] == 8:
            for mismatch in range(tuple[1]):
                ref_cursor += 1
                if ref_cursor not in tmp_variants:
                    tmp_variants[ref_cursor] = v.Variant(chromosome, ref_cursor)
                for segment in subsegments:
                    if segment.pos <= ref_cursor <= segment.end:
                        tmp_variants[ref_cursor].add_segment(segment.id, [seq[read_cursor], qual[read_cursor]])
                        break
                    if segment.end < ref_cursor:
                        del subsegments[subsegments.index(segment)]
                read_cursor += 1
        elif tuple[0] == 2:
            for deletion in range(tuple[1]):
                ref_cursor += 1
                if ref_cursor not in tmp_variants:
                    tmp_variants[ref_cursor] = v.Variant(chromosome, ref_cursor)
                for segment in subsegments:
                    if segment.pos <= ref_cursor <= segment.end:
                        tmp_variants[ref_cursor].add_segment(segment.id, ['-'])
                        break
                    if segment.end < ref_cursor:
                        del subsegments[subsegments.index(segment)]
        elif tuple[0] == 7:
            ref_cursor += tuple[1]
            read_cursor += tuple[1]
        elif tuple[0] == 1:
            read_cursor += tuple[1]


def remove_variations(start, end, chr):
    """
    Checks saved variants on position when all reads of that position are parsed. If a variant meets all selection 
    requirements it is kept, otherwise deleted
    :param start: 
    :param end: 
    :param chr: 
    """
    global tmp_variants, variants, segments_to_check
    for position in sorted(tmp_variants):
        cov = 0
        for segment in segments_to_check:
            if segments_to_check[segment][0] <= position <= segments_to_check[segment][1]:
                cov += 1
        if cov < NanoSV.opts_minimum_coverage:
            continue
        variant = tmp_variants[position]
        if position < start:
            continue
        elif position > end:
            break
        base_ratios = {'A': [0, 0], 'C': [0, 0], 'G': [0, 0], 'T': [0, 0], '=': [0, 0]}
        deletions = 0
        amount_of_vars = 0
        for segment_id, info in variant.segments.items():
            if info[0] == '-':
                deletions += 1
                amount_of_vars += 1
            elif info[1] >= NanoSV.opts_min_base_qual_ph:
                base_ratios[info[0]][0] += 1
                amount_of_vars += 1
            else:
                base_ratios[info[0]][1] += 1
                amount_of_vars += 1
        base_ratios['='][0] = cov - amount_of_vars
        qual_str = []
        if deletions < (NanoSV.opts_max_deletions * cov):
            haplotypes = sorted(base_ratios.items(), key=lambda x:sum(x[1]))[-2:]
            if haplotypes[0][0] != '=':
                for segment_id, info in variant.segments.items():
                    if info[0] == haplotypes[0][0]:
                        qual_str.append(str(info[1]))
            if haplotypes[1][0] != '=':
                for segment_id, info in variant.segments.items():
                    if info[0] == haplotypes[1][0]:
                        qual_str.append(str(info[1]))
            try:
                if haplotypes[0][1][0]/sum(haplotypes[0][1]) > NanoSV.opts_min_occurences_of_highq_var and haplotypes[1][1][0]/sum(haplotypes[1][1]) > NanoSV.opts_min_occurences_of_highq_var:
                    if sum(haplotypes[0][1])/(sum(haplotypes[1][1]) + sum(haplotypes[0][1])) > NanoSV.opts_min_occurences_of_var:
                        bin = int(position/NanoSV.opts_variant_bin_size)
                        variants[chr][bin][position] = tmp_variants[position]
            except ZeroDivisionError:
                ""
        del tmp_variants[position]


def write_bed():
    """
    Writes all phasing SNP positions to a bed output file.
    """
    with open("variants_of_sim.bed", 'w') as bedfile:
        for chromosome in variants:
            for bin in variants[chromosome]:
                for position, object in variants[chromosome][bin].items():
                    string = str(chromosome) + "\t" + str(position)
                    bedfile.write(string)