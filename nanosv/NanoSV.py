import configparser
import sys
import os
import argparse

import utils

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('bam', help='/path/to/file.bam')
parser.add_argument('-s','--sambamba', default='sambamba',type=str,help=' Give the full path to the sambamba or samtools executable [default: sambamba ]')
parser.add_argument('-c','--config', default=os.path.dirname(os.path.abspath(__file__))+"/config.ini",type=str,help='Give the full path to your own ini file [ default: config.ini ]')
parser.add_argument('-b','--bed', default=os.path.dirname(os.path.abspath(__file__))+"/bedfiles/human_hg19.bed",type=str,help=' Give the full path to your own bed file, used for coverage depth calculations [default: human_hg19.bed ]')
parser.add_argument('-o','--output',default=sys.stdout,type=argparse.FileType('w'),help='Give the full path to the output vcf file [default: <stdout> ]')
args = parser.parse_args()

opts_bam = args.bam    
opts_bed = args.bed
opts_sambamba = args.sambamba
opts_output = args.output    
cfg = configparser.ConfigParser()
if not args.config == os.path.dirname(os.path.abspath(__file__))+"/config.ini":
    cfg.read([ os.path.dirname(os.path.abspath(__file__))+"/config.ini", args.config ])
else:
    cfg.read( args.config )

opts_max_split = int(cfg.get('Filter options', 'max_split'))
opts_min_pid = float(cfg.get('Filter options', 'min_pid'))
opts_min_mapq = int(cfg.get('Filter options', 'min_mapq'))

opts_cluster_distance = int(cfg.get('Detection options', 'cluster_distance'))
opts_cluster_count = int(cfg.get('Detection options', 'cluster_count'))
opts_refread_distance = int(cfg.get('Detection options', 'refreads_distance'))
opts_hanging_length = int(cfg.get('Detection options', 'hanging_length'))
opts_mate_distance = int(cfg.get('Detection options', 'mate_distance'))
opts_depth_support = cfg.getboolean('Detection options', 'depth_support')

opts_qual_flag = int(cfg.get('Output filter options', 'qual_flag') )
opts_window_size = int(cfg.get('Output filter options', 'window_size'))
opts_svcluster = int(cfg.get('Output filter options', 'svcluster'))
opts_mapq_flag = int(cfg.get('Output filter options', 'mapq_flag'))
opts_pid_flag = float(cfg.get('Output filter options', 'pid_flag'))
opts_gap_flag = int(cfg.get('Output filter options', 'gap_flag'))
opts_ci_flag = int(cfg.get('Output filter options', 'ci_flag'))

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


    