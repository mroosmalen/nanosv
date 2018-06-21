NanoSV User Guide
=================

## Table of Contents
[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)
* [NanoSV](#nanosv)
  * [Summary](#summary)
  * [Installation](#installation)
  * [Citation](#citation)
* [Pre-processing](#pre-processing)
  * [Basecalling](#basecalling)
  * [Mapping](#mapping)
  * [LAST](#last-mapping)
    * [LAST installation](#last-installation)
    * [Running LAST](#running-last)
* [SV calling using NanoSV](#sv-calling-using-nanosv)
  * [NanoSV usage](#nanosv-usage)
  * [NanoSV arguments and parameters](#nanosv-arguments-and-parameters)
    * [Required arguments](#required-arguments)
    * [Optional arguments](#optional-arguments)
    * [Optional configuration parameters](#optional-configuration-parameters)
    * [Ancillary files that can be used for running NanoSV](#ancillary-files-that-can-be-used-for-running-nanosv)
  * [NanoSV output](#nanosv-output)
    
[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)

## NanoSV 
### Summary
NanoSV is a software package that can be used to identify structural genomic variations in long-read sequencing data, such as data produced by Oxford Nanopore Technologies’ MinION, GridION or PromethION instruments, or Pacific Biosciences RSII or Sequel sequencers.
NanoSV has been extensively tested using Oxford Nanopore MinION sequencing data, as described here: 
https://www.nature.com/articles/s41467-017-01343-4
https://link.springer.com/article/10.1007%2Fs00401-017-1743-5

The core algorithm of NanoSV identifies split- and gapped-aligned reads and clusters the reads according to the orientations and genomic positions of the read segments to define breakpoint-junctions of structural variations.

### Installation
NanoSV needs a working installation of python 3. You can install NanoSV using pip:
```
> pip install nanosv
```

Alternatively, you can get NanoSV from bioconda:
```
conda install -c bioconda nanosv 
```

### Citation
Cretu Stancu, M. *et al.* Mapping and phasing of structural variation in patient genomes using nanopore sequencing. Nat. Commun. 8, 1326 **(2017)**. (https://www.nature.com/articles/s41467-017-01343-4)

## Pre-processing

### Basecalling

Raw sequencing data can be basecalled using any available basecaller that is suited for your data. We use albacore for MinION/GridION/PromethION sequencing data, available through the nanopore community: http://community.nanoporetech.com/. Albacore can directly produce fastq files suitable for subsequent mapping to a reference genome.

### Mapping

NanoSV has been tested with different long read mappers, including BWA MEM, MINIMAP2, LAST and NGMLR. Below you can see how we tipically run these mappers. 

#### BWA MEM
```
> bwa mem -x ont2d -M -t 8 <reference> <fastq|fasta>
``` 
#### MINIMAP2
```
> minimap2 -t 8 -a <reference> <fastq|fasta>
``` 
#### NGMLR
```
> ngmlr -x ont -t 8 -r <reference> -q <fastq|fasta>
```
#### LAST 

We found that LAST alignments give the most accurate results for SV calling with NanoSV. However, mapping with LAST requires more compute resources. Follow the instructions below if you would like to use LAST alignment as input for your SV calling with NanoSV. 

##### LAST installation

Download the zip file from http://last.cbrc.jp/
```
> gunzip last.zip
> cd /path/to/lastdir
> make
```

##### Running LAST
First you need to index your reference genome by creating a lastal database:
```
> lastdb [referencedb] [reference.fa]
```
Train LAST to get the best scoring parameters for your particular alignment. We typically use a subset of > 10,000 reads for this step:
```
> last-train -Q1 [referencedb] [reads_sample.fastq] > [my-params]
```

Map your fastq data to reference:
```
> lastal -Q1 -p [my-params] [referencedb] [reads.fastq] | last-split > [reads.maf]
```
Convert the MAF file to SAM format:
```
> maf-convert -f [reference.dict] sam -r ‘ID:[id] PL:[nanopore] SM:[sample]’ [reads.maf] > [reads.sam]
```
The `[reference.dict]` file can be created by picard:
```
> java -jar picard.jar CreateSequenceDictionary REFERENCE=[reference.fa] OUTPUT=[reference.dict]
```
Convert SAM file to BAM file using sambamba (https://github.com/biod/sambamba) (samtools may function similarly):
```
> sambamba view -h -S --format=bam [reads.sam] > [reads.bam]
```
Sort the BAM file using sambamba: 
```
> sambamba sort [reads.bam]
```

All of the above commands can also be run at once using pipes:
```
> lastal -Q1 -p [my-params] [referencedb] [reads.fastq] | \
> last-split | \
> maf-convert -f [reference.dict] sam -r ‘ID:[id] PL:[nanopore] SM:[sample]’ /dev/stdin | \
> sambamba view -h -S --format=bam /dev/stdin | \
> sambamba sort /dev/stdin -o [reads.sorted.bam]

```

## SV calling using NanoSV

### NanoSV usage
```
> NanoSV [-h] [-t THREADS] [-s SAMBAMBA] [-c CONFIG] [-b BED] [-o OUTPUT] [reads.sorted.bam]
```

### NanoSV arguments and parameters:

#### required arguments:
```
bam              :   /path/to/reads.sorted.bam
```
This BAM file needs to be coordinate-sorted and indexed. Note that if you are performing SV calling on a large genome (e.g. human) and are only interested in calling intrachromosomal SVs, you may gain speed by splitting your BAM file by chromosome and running NanoSV per chromosome (on a compute cluster). 

#### optional arguments:
```
-h, --help       :   Show the help message and exit

-t, --threads    :   Maximum number of threads to use [default: 4 ]

-s, --sambamba   :   Give the full path to the sambamba or samtools executable [default: sambamba ]

-c, --config     :   Give the full path to your own ini file [ default: config.ini ]

-b, --bed        :   Give the full path to your own bed file, used for coverage depth calculations [default: human_hg19.bed ]

-o, --output     :   Give the full path to the output vcf file [default: <stdout> ]
```

#### optional configuration parameters:
NanoSV uses a config.ini file which contains default settings for all running parameters. Users can change the parameters by creating their own config.ini file and provide this as a command line argument [-c]
```
#Reads segments options
[Filter options]
# Maximum number of segments per read resulting from the mapping of the read to the reference sequence
max_split = 10
# Minimum percentage of identical bases of the mapped segment relative to the reference sequence
min_pid = 0.7
# Minimum mapping quality of the segment
min_mapq = 20

#Parameters for tuning detection and clustering of breakpoints:
[Detection options]
# Maximum distance between two adjacent break-end positions
cluster_distance = 10
# Minimum number of breakpoint-junctions (i.e. split-read junctions) for clustering
cluster_count = 2
# Minimum flanking length, to consider a read a reference read
refreads_distance = 100
# Minimum length of unmapped sequence for hanging reads that overlap a break-end
hanging_length = 20
# Maximum distance to search for the MATEID, i.e. a reciprocal breakpoint-junction, for example an inversion consist of two breakpoint-junctions (3’-to-3’ and 5’-to-5’)
mate_distance = 300
# If True, NanoSV will check the depth of coverage for possible breakpoint-junctions with orientations that indicate a possible deletion or duplication (3’-to-5’ and 5’-to-3’)
depth_support = True
# Minimum indel size to call gap and create subsegments
min_indel_size = 30

#Parameters for setting the FILTER flag in the vcf output:
[Output filter options]
# Filter flag: LowQual, if the QUAL score is lower
qual_flag = 20
# Filter flag: SVcluster, if there are more SVs within a window size, they will be marked as SVcluster
window_size = 1000
# Filter flag: SVcluster, indicating the number of SVs within a certain window size (set by window_size above)
svcluster = 2
# Filter flag: MapQual, if the median mapq is lower than specified by this parameter
mapq_flag = 80
# Filter flag: PID, if the median percentage identity is lower than specified by this parameter
pid_flag = 0.80
# Filter flag: Gap, if the median GAP is higher than specified by this parameter
gap_flag = 100
# Filter flag: CIPOS|CIEND, if the CIPOS|CIEND bigger than specified by this parameter
ci_flag = 30

[Phasing Options]
##Phasing is still experimental for now and needs to be properly tested and benchmarked, so use it under your own responsibility. 
#If True, NanoSV will use phasing as an addition in calling SVs
phasing_on = False
#SNP positions are stored in bins to improve speed. This setting sets the bin size
variant_bin_size = 1000000
#Window measured from the breakpoint in which SNPs are sought to be used in read clustering
phasing_window = 7000
#Minimum coverage to call a SNP for phasing
min_coverage = 10
#Maximum percentage of deletions on position
max_deletions = 0.25
#Minimum occurence of variant to call a SNP for phasing
min_occurences_of_var = 0.4
#Minimum occurence of high quality calls of certain variant
min_highq_var = 0.6
#Minimum quality to call 'high quality'
min_base_qual_ph = 11
#cut-off setting to stop clustering if highest similarity between reads is too low
clustering_cutoff = 0.3
```

#### Ancillary files that can be used for running NanoSV:
To estimate a coverage increase or decrease near predicted breakpoint-junctions, the average coverage across a putative deletion or duplication interval is compared to the distribution of coverage across random positions in the reference sequence. This calculation is only performed if `depth_support = True` in config.ini. A default bed file is provided that contains 1,000,000 random positions on the hg19/GRCh37 human genome reference, excluding simple repeat regions (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz) and gap regions (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz). The file format is standard BED format (chr\<TAB\>startpos\<TAB\>endpos).

#### NanoSV output
NanoSV will output a VCF file. After dealing with mis-sinterpration of the breakpoint orientation in the past, we decided to call only breakpoints instead of SV types (such as inversions, deletions, etc.) We want to leave the interpretation of the breakpoints to our users and leave the output as "assumption-free" as possible.
These breakpoints (BNDs) are reported using the standard VCF specifications described in https://samtools.github.io/hts-specs/VCFv4.2.pdf (chapter 5.4) .
Using the depth_support mode (on by default), we test if there is a significant coverage change around a breakpoint with the right orientation to be able to interpret deletions and duplications.
The reported breakpoints will also have flags in the FILTER field according to the threshold filter values chosen in the config file.

### Phasing on NanoSV
We are working on haplotype-aware SV detection with NanoSV. However, it still needs to be properly tested and benchmarked. If you feel brave enough to try it, go for it (and report back to us with feedback). It will try to look for useful phasing SNPs around the breakpoints on the long read data, increasing a lot the runtime. Please check the appropiate settings on the config file.
