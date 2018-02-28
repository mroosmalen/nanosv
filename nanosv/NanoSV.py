import configparser
import sys
import os
import argparse

import utils

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('bam', help='/path/to/file.bam')
parser.add_argument('-s','--sambamba', default='sambamba',type=str,help='Path to sambamba')
parser.add_argument('-c','--config', default=os.path.dirname(os.path.abspath(__file__))+"/config.ini",type=str,help='Path to config.ini file')
parser.add_argument('-b','--bed', default=os.path.dirname(os.path.abspath(__file__))+"/bedfiles/hg19_1M_random_protein_coding.bed",type=str,help='Path to bed file')
parser.add_argument('-o','--output',default=sys.stdout,type=argparse.FileType('w'),help='Path to output file')
args = parser.parse_args()

opts_bam = args.bam    
opts_bed = args.bed
opts_sambamba = args.sambamba
opts_output = args.output    
cfg = configparser.ConfigParser()
cfg.read(args.config)

opts_max_split = int(cfg.get('Filter Options', 'max_split'))
opts_min_pid = float(cfg.get('Filter Options', 'min_pid'))
opts_min_mapq = int(cfg.get('Filter Options', 'min_mapq'))

opts_cluster_distance = int(cfg.get('Detection Options', 'cluster_distance'))
opts_cluster_count = int(cfg.get('Detection Options', 'cluster_count'))
opts_refread_distance = int(cfg.get('Detection Options', 'refread_distance'))
opts_hanging_length = int(cfg.get('Detection Options', 'hanging_length'))
opts_mate_distance = int(cfg.get('Detection Options', 'mate_distance'))
opts_depth_support = cfg.getboolean('Detection Options', 'depth_support')

opts_qual_flag = int(cfg.get('Output Filter Options', 'qual_flag') )
opts_window_size = int(cfg.get('Output Filter Options', 'window_size'))
opts_svcluster = int(cfg.get('Output Filter Options', 'svcluster'))
opts_mapq_flag = int(cfg.get('Output Filter Options', 'mapq_flag'))
opts_pid_flag = float(cfg.get('Output Filter Options', 'pid_flag'))
opts_gap_flag = int(cfg.get('Output Filter Options', 'gap_flag'))
opts_ci_flag = int(cfg.get('Output Filter Options', 'ci_flag'))

def main():    
    if opts_depth_support:
        utils.coverage.calculate_coverage_bed()
    else:
        coverages = []
    
    utils.parse_bam.parse_bam()
    
    utils.create_vcf.print_vcf_header()
    
    utils.parse_reads.parse_reads()
    
    utils.parse_breakpoints.parse_breakpoints()

    utils.parse_svs.parse_svs()    

if __name__ == "__main__":
    main()


    