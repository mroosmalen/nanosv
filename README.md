# nanosv

Nanosv is a structural variation detection tool for Oxford Nanopore data. Nanosv works for now only bam files, which are mapped with LAST.

# Download/Install
```
git clone --recursive https://github.com/mroosmalen/nanosv
```

# How to run
```
perl nanosv.pl [ options ] sample.bam > svs.vcf
```

# Options
```
  General options:
  -h|help                       Help
  -t|threads            [i]     Number of threads. Default: 8
  -sambamba             [s]     Path to sambamba: Default: sambamba_v0.6.3
  
  Filter options:
  -s|split              [i]     Maximum number of segments per read. Default: 10
  -p|pid                [s]     Minimum percentage identity to reference. Default: 0.70
  -m|mapq               [i]     Minimum mapping qualty. Default: 20  
  
  Detection options:
  -d|distance           [i]     Maximum distance to cluster SVs together. Default: 10  
  -c|count              [i]     Minimum number of supporting reads. Default: 2  
  -f|refdistance        [i]     Minimum distance for reference reads: Default: 100
  -u|unmapped           [i]     Minimum unmapped length of hanging segments: 20
  -r|matedistance       [i]     Maximum distance to look for mateid. Default: 300  
    
  Output filter options:
  -w|window             [i]     Maximum window size. Default: 1000
  -n|cluster            [i]     Maximum number of SV's in a window. Default: 2
  -q|mapqf              [i]     Minimum median mapping quality of a SV. Default: 80
  -i|pidf               [s]     Minimum median percentage identity to reference. Default: 0.80
  -g|gap                [i]     Maximum median gap size. Default: 100  
  -y|ci                 [i]     Maximum Confidence interval distance. Default: 20
```
