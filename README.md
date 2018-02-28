NanoSV User Guide
=================

## Table of Contents
[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)
* [NanoSV summary](#nanosv-summary)
* [Pre-processing](#pre-processing)
  * [Basecalling](#basecalling)
  * [LAST mapping](#last-mapping)
    * [LAST installation](#last-installation)
    * [Running LAST](#running-last)
* [SV calling using NanoSV](#sv-calling-using-nanosv)
  * [NanoSV installation](#nanosv-installation)
  * [NanoSV usage](#nanosv-usage)
  * [NanoSV arguments and parameters](#nanosv-arguments-and-parameters)
    * [Required arguments](#required-arguments)
    * [Optional arguments](#optional-arguments)
    * [Optional configuration parameters](#optional-configuration-parameters)
    
[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)

## NanoSV summary
NanoSV is a software package that can be used to identify structural genomic variations in long-read sequencing data, such as data produced by Oxford Nanopore Technologies’ MinION, GridION or PromethION instruments, or Pacific Biosciences sequencers.
NanoSV has been extensively tested using Oxford Nanopore MinION sequencing data, as described here: https://www.nature.com/articles/s41467-017-01343-4
The core algorithm of NanoSV identifies split mapped reads and clusters the split-mapped orientations and genomic positions to identify breakpoint-junctions of structural variations.

## Pre-processing

### Basecalling

Raw sequencing data can be basecalled using any available basecaller that is suited for your data. We use albacore for MinION/GridION/PromethION sequencing data, available through the nanopore community: http://community.nanoporetech.com/. Albacore can directly produce fastq files suitable for subsequent mapping to a reference genome.

### LAST mapping

The current version of NanoSV uses LAST mapping data. Newer releases are likely also compatible with alignments of other mappers, e.g. BWA MEM.

#### LAST installation

Download the zip file from http://last.cbrc.jp/
```
> gunzip last.zip
> cd /path/to/lastdir
> make
```

#### Running LAST
First you need to index your reference genome by creating a lastal database:
```
> lastdb [referencedb] [reference.fa]
```
Train LAST to get the best scoring parameters for the alignment. We typically use > 10,000 reads for this step
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
Convert SAM file to BAM file using sambamba (samtools may function similarly)
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

### NanoSV installation
```
> pip install nanosv
```

### NanoSV usage
```
> NanoSV [-h] [-s SAMBAMBA] [-c CONFIG] [-b BED] [-o OUTPUT] [reads.sorted.bam]
```

### NanoSV arguments and parameters:

#### required arguments:
bam /path/to/reads.sorted.bam

#### optional arguments:
-h, --help       :   show this help message and exit

-s, --sambamba   :   Give the full path to sambamba executable [default: sambamba ]

-c, --config     :   Give the full path to the ini file [ default: config.ini ]

-b, --bed BED    :   Give the full path to the bed file [default: default.bed ]

-o, --output     :   Give the full path to the output vcf file [default: stdout ]


#### optional configuration parameters:
NanoSV uses a config.ini file which contains default settings for all running parameters. Users can change the parameters by creating their own config.ini file and provide this as a command line argument [-c]
```
#Reads segments options
[Filter options]
# Maximum number of segments per read resulting from the mapping of the read the a reference sequence
max_split = 10
# Minimum percentage of identical bases of the mapped segment relative to the reference sequence      
min_pid = 0.7
# Minimum mapping quality of the segment
min_mapq = 20
```
