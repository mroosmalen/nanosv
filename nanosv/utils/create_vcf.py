#!/usr/bin/python

import io
import sys
import time
import vcf as py_vcf
import os

from utils import parse_bam as bam
from version import __version__

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NanoSV

def print_vcf_header():
    """
    Creates vcf header by setting vcf format and calling functions to create each header section.
    """
    global vcf_writer
    
    vcf = py_vcf.Reader(io.StringIO("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+ str(bam.sample_name)))
    set_metadata_header(vcf)
    set_contig_header(vcf)
    set_info_header(vcf)
    set_alt_header(vcf)
    set_format_header(vcf)
    set_filter_header(vcf)
    vcf_writer = py_vcf.Writer(NanoSV.opts_output, vcf)

def set_metadata_header(vcf):
    """
    Create metadata header in vcf output file
    :param vcf used to add contigs to file:
    """
    vcf.metadata = { 'fileformat': 'VCFv4.1',
                     'fileDate': str(time.strftime("%c")),
                     'cmdline': [" ".join(sys.argv)]
                   }

def set_contig_header(vcf):
    """
    Create contig header in vcf output file
    :param vcf used to add contigs to file:
    """
    for item in bam.header["SQ"]:
        vcf.contigs[ item['SN'] ] = py_vcf.parser._Contig(item["SN"], item["LN"])

def set_info_header(vcf):
    """
    Create info header in vcf output file
    :param vcf used to add infos to file:
    """
    
    vcf.infos = { 'IMPRECISE': py_vcf.parser._Info("IMPRECISE", 0, "Flag", "Imprecise structural variant", "NanoSV", __version__),
                  'SVTYPE': py_vcf.parser._Info("SVTYPE", 1, "String", "Type of structural variant", "NanoSV", __version__),
                  'SVMETHOD': py_vcf.parser._Info("SVMETHOD", 1, "String", "Type of approach used to detect SV", "NanoSV", __version__),
                  'END': py_vcf.parser._Info("END", 1, "Integer", "End position of structural variant", "NanoSV", __version__),
                  'CIPOS': py_vcf.parser._Info("CIPOS", 2, "Integer", "Confidence interval around POS", "NanoSV", __version__),
                  'CIEND': py_vcf.parser._Info("CIEND", 2, "Integer", "Confidence interval around END", "NanoSV", __version__),
                  'SVLEN': py_vcf.parser._Info("SVLEN", None, "Integer", "Distance between the two genomic positions", "NanoSV", __version__),
                  'RT': py_vcf.parser._Info("RT", 3, "Integer", "Number of the different read types (2d, template, complement)","NanoSV", __version__),
                  'GAP': py_vcf.parser._Info("GAP", 1, "Integer","Median number of bases between the two segments of the SV, in case of an insertion this is the size of the insertion","NanoSV", __version__),
                  'MAPQ': py_vcf.parser._Info("MAPQ", 2, "Integer","Median mapping quality of the two segments of the structural variant", "NanoSV",__version__),
                  'PID': py_vcf.parser._Info("PID", 2, "Float","Median percentage identity to the reference of the two segments of the structural variant","NanoSV", __version__),
                  'PLENGTH': py_vcf.parser._Info("PLENGTH", 2, "Float","Median segment length percentage of the two segments of the structural variant","NanoSV", __version__),
                  'RLENGTH': py_vcf.parser._Info("RLENGTH", 1, "Integer", "Median length of the total reads", "NanoSV", __version__),
                  'MATEID': py_vcf.parser._Info('MATEID', None, 'String', 'ID of mate breakend', 'NanoSV', __version__)
                }
    if NanoSV.opts_depth_support:
        vcf.infos['DEPTHPVAL'] = py_vcf.parser._Info("DEPTHPVAL", 1, "Float","In case of a duplication or deletion the P-value of the significance test is shown here","NanoSV", __version__)    

def set_format_header(vcf):
    """
    Create format header in vcf output file
    :param vcf used to add formats to file:
    """
    vcf.formats = { 'GT': py_vcf.parser._Format('GT', 1, 'String', 'Genotype'),
                    'DR': py_vcf.parser._Format('DR', 2, 'Integer', 'Number of reference reads'),
                    'DV': py_vcf.parser._Format('DV', 2, 'Integer', 'Number of variant reads'),
                    'GQ': py_vcf.parser._Format('GQ', 1, 'Integer', 'Genotype quality'),
                    'HR': py_vcf.parser._Format('HR', 2, 'Integer', 'Number of hanging variant reads'),
                    'PL': py_vcf.parser._Format('PL', 'G', 'Integer','Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification')
                  }

def set_filter_header(vcf):
    """
    Create filter header in vcf output file
    :param vcf used to add filters to file:
    """
    vcf.filters = { 'SVcluster': py_vcf.parser._Filter("SVcluster","There are more than " + str(NanoSV.opts_svcluster) + " SVs in a window of " + str(NanoSV.opts_window_size) + " on both sides"),
                    'Gap': py_vcf.parser._Filter("Gap", "The median gap size is larger than " + str(NanoSV.opts_gap_flag) + " for non insertions"),
                    'MapQual': py_vcf.parser._Filter("MapQual", "The median mapping quality is less than " + str(NanoSV.opts_mapq_flag)),
                    'PID': py_vcf.parser._Filter("PID", "The median PID of one of the two sides is less than " + str(NanoSV.opts_pid_flag)),
                    'CIPOS': py_vcf.parser._Filter("CIPOS", "The CIPOS is greater than " + str(NanoSV.opts_ci_flag)),
                    'CIEND': py_vcf.parser._Filter("CIEND", "The CIEND is greater than " + str(NanoSV.opts_ci_flag)) ,
                    'LowQual': py_vcf.parser._Filter('LowQual','QUAL score is less than 20')
                  }

def set_alt_header(vcf):
    """
    Create alt header in vcf output file.
    :param vcf used to add alts to file:
    """
    vcf.alts = { 'DEL': py_vcf.parser._Alt("DEL", "Deletion"), 
                 'DUP': py_vcf.parser._Alt("DUP", "Duplication"),
                 'BND': py_vcf.parser._Alt("BND", "Breakend"), 
                 'INS': py_vcf.parser._Alt("INS", "Insertion")
               }
