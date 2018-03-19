#!/usr/bin/python
import sys

from math import log10
from statistics import median
from version import __version__

import collections
import math
import os
import vcf as py_vcf

from utils import parse_bam as bam
from utils import coverage

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

class SV:
    def __init__(self, id, breakpoint):
        self.id = id
        self.chr = breakpoint.segment_1["rname"]
        self.chr2 = breakpoint.segment_2["rname"]
        self.flag1 = breakpoint.segment_1["flag"]
        self.flag2 = breakpoint.segment_2["flag"]
        self.pos = [breakpoint.segment_1["pos"]]
        self.ref = 'N'
        self.alt = str(breakpoint.svtype)
        self.qual = None
        self.filter = []
        self.info = {
            'IMPRECISE': True,
            'END': [breakpoint.segment_2["pos"]],
            'SVTYPE': breakpoint.svtype,
            'SVMETHOD': 'NanoSV_v'+__version__,
        }
        self.format = {
            'GT': './.',
            'DV': [1, 1],
            'VO': [(1 - 10 ** (-breakpoint.segment_1["mapq"] / 10.0)),
                   (1 - 10 ** (-breakpoint.segment_2["mapq"] / 10.0))],
            'DR': [0, 0],
            'RO': [0, 0],
            'HR': [0, 0]
        }
        self.breakpoints = [breakpoint.id]
        self.set = 0
        self.SVcluster = 1

    def addBreakpoint(self, breakpoint):
        """
        Adds breakpoint to SV object
        :param breakpoint is breakpoint-object:
        """
        self.breakpoints.append(breakpoint.id)
        self.pos.append(breakpoint.segment_1["pos"])
        self.info['END'].append(breakpoint.segment_2["pos"])
        self.format['DV'][0] += 1
        self.format['DV'][1] += 1
        self.format['VO'][0] += (1 - 10 ** (-breakpoint.segment_1["mapq"] / 10.0))
        self.format['VO'][1] += (1 - 10 ** (-breakpoint.segment_2["mapq"] / 10.0))

    def addInfoField(self, key, value):
        """
        Adds info field to SV object
        :param key is the name of the info field result:
        :param value is the value of the field result:
        """
        if key in self.info:
            if isinstance(self.info[key], list):
                try:
                    if isinstance(self.info[key][0], list) or isinstance(self.info[key][1], list):
                        if value[0] is not None: self.info[key][0].append(value[0][0])
                        if value[1] is not None: self.info[key][1].append(value[1][0])
                    else:
                        self.info[key].append(value[0])

                except IndexError:
                    self.info[key].append(value[0])
        else:
            self.info[key] = value

    def significanceTest(self, avg_dupdel_cov, dup):
        """
        Calculates if there is a significant difference in coverage to call DUP or DEL
        :param avg_dupdel_cov is the distribution with coverage of random positions in the genome:
        :param dup is a boolean. True in case of a duplication, False in case of a deletion:
        :return True or false depending on the outcome of the test. True when there is a significant difference:
        """
        teller = 0
        if dup:
            for value in coverage.coverages:
                if float(value) > avg_dupdel_cov:
                    teller += 1
        else:
            for value in coverage.coverages:
                if float(value) < avg_dupdel_cov:
                    teller += 1
        if teller / len(coverage.coverages) < 0.05:
            self.info["DEPTHPVAL"] = '%.3f' % (teller / len(coverage.coverages))
            return True
        else:
            self.info["DEPTHPVAL"] = '%.3f' % (teller / len(coverage.coverages))
            return False
    
    def setReferenceBase(self):
        """
        Get the reference bases for the SV object
        :param
        """
        bases = dict()
        for pileupcolumn in bam.bam.pileup(self.chr,self.pos-1,self.pos, truncate=True):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    if not pileupread.alignment.query_sequence[pileupread.query_position] in bases:
                        bases[ pileupread.alignment.query_sequence[pileupread.query_position] ] = 0
                    bases[ pileupread.alignment.query_sequence[pileupread.query_position] ] += 1
        if ( len(bases.keys()) > 0 ):
            self.ref = sorted(bases, key=bases.__getitem__)[-1]   
    
    def setArguments(self, opts_depth_support):
        """
        Set standard info field arguments.
        :param opts_depth_support is a boolean for if the coverage dup/del analysis is switched on:
        """
        self.info['CIPOS'] = [ int(min(self.pos) - median(self.pos)), int(max(self.pos) - median(self.pos)) ]
        self.info['CIEND'] = [ int(min(self.info['END']) - median(self.info['END'])), int(max(self.info['END']) - median(self.info['END'])) ]
        if self.info['CIPOS'] == [0,0] and self.info['CIEND'] == [0,0]:
            del self.info['IMPRECISE']
        self.pos = int(median(self.pos))
        self.info['END'] = int(median(self.info['END']))
        if self.chr == self.chr2:
            self.info['SVLEN'] = (self.info['END'] - self.pos)
        self.setInfoField()
        self.setReferenceBase()
        dup = 0
        if self.alt == "INS":
            self.alt = py_vcf.model._SV("INS")
        if self.alt == "BND":
            if self.flag1 & 16:
                if self.flag2 & 16:
                    self.alt = py_vcf.model._Breakend(self.chr2, int(self.info['END']), True, False, self.ref, True)
                    if opts_depth_support:
                        avg_dupdel_cov = self.getDupDelcoverage()
                        if self.significanceTest(avg_dupdel_cov, True):
                            self.alt = py_vcf.model._SV("DUP")
                            self.info['SVTYPE'] = "DUP"
                        dup = 1
                else:
                    self.alt = py_vcf.model._Breakend(self.chr2, int(self.info['END']), True, True, self.ref, True)
            else:
                if self.flag2 & 16:
                    self.alt = py_vcf.model._Breakend(self.chr2, int(self.info['END']), False, False, self.ref, True)
                else:
                    self.alt = py_vcf.model._Breakend(self.chr2, int(self.info['END']), False, True, self.ref, True)
                    if opts_depth_support:
                        avg_dupdel_cov = self.getDupDelcoverage()
                        if self.significanceTest(avg_dupdel_cov, False):
                            self.alt = py_vcf.model._SV("DEL")
                            self.info['SVTYPE'] = "DEL"
        gt_lplist = self.bayes_gt(sum(self.format['RO']), sum(self.format['VO']), dup)
        gt_sum = 0
        
        pplist = []
        for gt in gt_lplist:
            pplist.append( 10 ** gt )
            gt_sum += 10 ** gt
        
        if gt_sum > 0:
            npplist = [x / gt_sum for x in pplist]
            
            plplist = []
            for p in npplist:
                if p > 0:
                    plplist.append( int( -10 * log10( p ) ) )
                else:
                    plplist.append( 9999 )
            
            gt_idx = plplist.index(min(plplist))
            gq_idx = plplist.index(sorted(plplist)[1])
            
            plplist[gt_idx] = 0
            
            if gt_idx == 0:
                self.format['GT'] = '0/0'
            elif gt_idx == 1:
                self.format['GT'] = '0/1'
            elif gt_idx == 2:
                self.format['GT'] = '1/1'
            
            self.format['PL'] = plplist
            self.format['GQ'] = plplist[gq_idx]
            self.qual = plplist[0]
        else:
            self.format['GT'] = './.'
            self.format['PL'] = [0, 0, 0]
            self.format['GQ'] = [0]
            self.qual = 0
            
        self.set = 1

    def getDupDelcoverage(self):
        """
        Returns the average coverage of all positions within the start and end of the SV.
        :return average coverage:
        """
        start = round(self.pos)
        stop = round(self.info['END'])
        if start > stop:
            start = round(self.info['END'])
            stop = round(self.pos)
        position = str(self.chr) + ":" + str(start) + "-" + str(stop)
        dupdel_coverages = []
        
        if 'sambamba' in NanoSV.opts_sambamba:
            cmd = NanoSV.opts_sambamba + " depth base --min-coverage=0 " + NanoSV.opts_bam + " -L " + position + " 2> /dev/null | awk '{if (NR!=1) print $3}'"
        elif 'samtools' in NanoSV.opts_sambamba:
            cmd = NanoSV.opts_sambamba + " depth " + NanoSV.opts_bam + " -r " + position + " | awk '{if (NR!=1) print $3}'"
        
        with os.popen(cmd) as commandoutput:
            for line in commandoutput:
                dupdel_coverages.append(float(line.rstrip()))
        return float(sum(dupdel_coverages)) / len(dupdel_coverages)

    def log_choose(self, n, k):
        """
        Returns log of given variable.
        :param n:
        :param k:
        """
        r = 0.0
        # swap for efficiency if k is more than half of n
        if k * 2 > n:
            k = n - k
        for d in range(1, k + 1):
            r += math.log(n, 10)
            r -= math.log(d, 10)
            n -= 1

        return (r)

    def bayes_gt(self, ref, alt, is_dup):
        """
        Returns genotype list
        :param ref:
        :param alt:
        :param is_dup boolean for if SV is dup:
        :return list with probabilities of different genotypes:
        """
        # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
        if is_dup:  # specialized logic to handle non-destructive events such as duplications
            #p_alt = [1e-2, 1 / 3.0, 0.5]
            p_alt = [1e-3, 0.5, 0.9]
        else:
            p_alt = [1e-3, 0.5, 0.9]

        total = ref + alt
        log_combo = self.log_choose(int(total), int(alt))

        lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
        lp_het = log_combo + (alt * math.log(p_alt[1], 10)) + (ref * math.log(1 - p_alt[1], 10))
        lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

        return [lp_homref, lp_het, lp_homalt]

    def setInfoField(self):
        """
        Calculates RT, MAPQ, PID and PLENGTH
        """
        for field in self.info:
            if field == "RT":
                rt = [0, 0, 0]
                for type in self.info[field]:
                    if type == "2d": rt[0] += 1
                    if type == "template": rt[1] += 1
                    if type == "complement": rt[2] += 1
                self.info[field] = rt
            elif isinstance(self.info[field], list):
                if isinstance(self.info[field][0], list):
                    self.info[field][0] = median(map(float, self.info[field][0]))
                    self.info[field][1] = median(map(float, self.info[field][1]))
                    if field == "MAPQ":
                        self.info[field][0] = int(self.info[field][0])
                        self.info[field][1] = int(self.info[field][1])
                    elif field == "PID" or field == "PLENGTH":
                        self.info[field][0] = round(self.info[field][0], 3)
                        self.info[field][1] = round(self.info[field][1], 3)
                else:
                    if not 'CI' in field:
                        self.info[field] = self.median_type(self.info[field])

    def median_type(self, toCalculate):
        """
        Returns median of given list in the original datatype
        :param toCalculate list to calculate median for:
        :return median
        """
        if type(toCalculate[0]) == float:
            return median(toCalculate)
        elif type(toCalculate[0]) == int:
            return int(median(toCalculate))

    def setFormat(self):
        """
        Creates format line of VCF result
        :return a list with the order of the values and the values:
        """
        del self.format['VO']
        del self.format['RO']
        format_list = sorted(self.format.keys())
        format_list.remove('GT')
        format_list.insert(0, 'GT')
        format = ':'.join(format_list)
       
        call = py_vcf.model._Call('site', bam.sample_name, collections.namedtuple('CallData', format_list)(**self.format))
        samples_indexes = [0]
        samples = [call]
        return [format, samples_indexes, samples]

    def printVCF(self, vcf_writer):
        """
        Combines all the values in a pyvcf record to write to the vcf
        :param vcf_writer is necessary to write a record to the vcf:
        """
        format_output = self.setFormat()
        if isinstance(self.alt, py_vcf.model._Breakend):
            del self.info['END']
        record = py_vcf.model._Record(self.chr, self.pos, self.id, self.ref, [self.alt], self.qual, self.filter,
                                      self.info, format_output[0], format_output[1], format_output[2])
        vcf_writer.write_record(record)
