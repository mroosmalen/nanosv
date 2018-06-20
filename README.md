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
  * [LAST mapping](#last-mapping)
    * [LAST installation](#last-installation)
    * [Running LAST](#running-last)
* [SV calling using NanoSV](#sv-calling-using-nanosv)
  * [NanoSV usage](#nanosv-usage)
  * [NanoSV arguments and parameters](#nanosv-arguments-and-parameters)
    * [Required arguments](#required-arguments)
    * [Optional arguments](#optional-arguments)
    * [Optional configuration parameters](#optional-configuration-parameters)
    * [Ancillary files that can be used for running NanoSV](#ancillary-files-that-can-be-used-for-running-nanosv)
    
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

NanoSV has been tested with different long read mappers, including BWA MEM, MINIMAP2, LAST and NGMLR.

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
#### LAST mapping

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
> NanoSV [-h] [-s SAMBAMBA] [-c CONFIG] [-b BED] [-o OUTPUT] [reads.sorted.bam]
```

### NanoSV arguments and parameters:

#### required arguments:
```
bam              :   /path/to/reads.sorted.bam
```
Note that if you are performing SV calling on a large genome (e.g. human) and are only interested in calling intrachromosomal SVs, you may gain speed by splitting your BAM file by chromosome and running NanoSV per chromosome (on a compute cluster).

#### optional arguments:
```
-h, --help       :   Show the help message and exit

-s, --sambamba   :   Give the full path to the sambamba or samtools executable [default: sambamba ]

-c, --config     :   Give the full path to your own ini file [ default: config.ini ]

-b, --bed        :   Give the full path to your own bed file, used for coverage depth calculations [default: human_hg19.bed ]

-o, --output     :   Give the full path to the output vcf file [default: <stdout> ]
```

#### optional configuration parameters:
NanoSV uses a config.ini file which contains default settings for all running parameters. Users can change the parameters by creating their own config.ini file and provide this as a command line argument [-c]
```
#Reads and segments options
[Filter options]
# Maximum number of segments per read resulting from the mapping of the read the a reference sequence
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
# Minimum flanking sequence length, to consider a read a reference read, i.e. the sequence mapped left and right of the breakpoint should be larger than the set value
refreads_distance = 100
# Minimum length of unmapped sequence for including reads as hanging reads that overlap (support) a break-end
hanging_length = 20
# Maximum distance to search for the MATEID, i.e. a reciprocal breakpoint-junction, for example an inversion consist of two breakpoint-junctions (3’-to-3’ and 5’-to-5’)
mate_distance = 300
# If TRUE, NanoSV will check the depth of coverage for possible breakpoint-junctions with orientations that indicate a possible deletion or duplication (3’-to-5’ and 5’-to-3’). Needs an auxiliar bed file, provided with -b to the main NanoSV command.
depth_support = True

#Parameters for setting the FILTER flag in the vcf output:
[Output filter options]
# Filter flag: LowQual, set if the QUAL score of the called structural variation is lower
qual_flag = 20
# Filter flag: SVcluster, set if there are more SVs within a window size, they will be marked as SVcluster
window_size = 1000
# Filter flag: SVcluster, set if the number of SVs within a certain window size (set by window_size above) exceeds this treshold
svcluster = 2
# Filter flag: MapQual, set if the median mapq is lower than specified by this parameter
mapq_flag = 80
# Filter flag: PID, set if the median percentage identity is lower than specified by this parameter
pid_flag = 0.80
# Filter flag: Gap, set if the median GAP is higher than specified by this parameter
gap_flag = 100
# Filter flag: CIPOS|CIEND, set if the CIPOS|CIEND is larger than specified by this parameter
ci_flag = 30

[Runtime options]
# Number of threads to use with sambamba
sambamba_threads = 1
```

#### Ancillary files that can be used for running NanoSV:
To estimate a coverage increase or decrease near predicted breakpoint-junctions, the average coverage across a putative deletion or duplication interval is compared to the distribution of coverage across random positions in the reference sequence. This calculation is only performed if `depth_support = True` in config.ini. A default bed file is provided that contains 1,000,000 random positions on the hg19/GRCh37 human genome reference, excluding simple repeat regions (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz) and gap regions (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz). The file format is standard BED format (chr\<TAB\>startpos\<TAB\>endpos).

