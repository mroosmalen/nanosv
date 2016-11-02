#! /usr/bin/perl

use strict;
use Getopt::Long;
use List::Util qw( min max sum );
use Data::Dumper;
use FindBin;
use POSIX qw(strftime);

use lib "$FindBin::Bin";
use Read;
use Segment;
use Breakpoint;
use SV;

sub usage {
  my ( $message ) = @_;
  
  warn $message . "\n" if $message;
  
  warn <<END;
  
  Usage:
   
  Run by typing:	perl nanosv.pl [options] <.bam>
  
  General options:
  -h|help			Help
  -t|threads		[i]	Number of threads. Default: 8
  -sambamba		[s]	Path to sambamba: Default: sambamba_v0.6.3
  
  Filter options:
  -s|split		[i]	Maximum number of segments per read. Default: 10
  -p|pid		[s]	Minimum percentage identity to reference. Default: 0.70
  -m|mapq		[i]	Minimum mapping qualty. Default: 20  
  
  Detection options:
  -d|distance		[i]	Maximum distance to cluster SVs together. Default: 10  
  -c|count		[i]	Minimum number of supporting reads. Default: 2  
  -f|refdistance	[i]	Minimum distance for reference reads: Default: 100
  -u|unmapped 		[i]	Minimum unmapped length of hanging segments: 20
  -r|matedistance	[i]	Maximum distance to look for mateid. Default: 300  
    
  Output filter options:
  -w|window		[i]	Maximum window size. Default: 1000
  -n|cluster		[i]	Maximum number of SV's in a window. Default: 2
  -q|mapqf		[i]	Minimum median mapping quality of a SV. Default: 80
  -i|pidf		[s]	Minimum median percentage identity to reference. Default: 0.80
  -g|gap		[i]	Maximum median gap size. Default: 100  
  -y|ci			[i]	Maximum Confidence interval distance. Default: 20
  
END
  exit;
}

my %opt;
%opt = (
  'help'		=> undef,
  'bam'			=> undef,
  'threads'		=> 8,
  'split'		=> 10,
  'distance'		=> 10,
  'matedistance'	=> 300,
  'window'		=> 1000,
  'refdistance'		=> 100,
  'unmapped'		=> 20,
  'cluster'		=> 2,
  'count'		=> 2,
  'pid'			=> 0.70,
  'mapq'		=> 20,
  'mapqf'		=> 80,
  'pidf'		=> 0.80,
  'gap'			=> 100,
  'ci'			=> 20,
  'sambamba'		=> 'sambamba_v0.6.3'
);

GetOptions (
  'h|help'		=> \$opt{help},
  't|threads=i'		=> \$opt{threads},
  's|split=i'		=> \$opt{split},
  'd|distance=i'	=> \$opt{distance},
  'r|matedistance=i'	=> \$opt{matedistance},
  'f|refdistance=i'	=> \$opt{refdistance},
  'u|unmapped=i'	=> \$opt{unmapped},
  'w|window=i'		=> \$opt{window},
  'n|cluster=i'		=> \$opt{cluster},
  'c|count=i'		=> \$opt{count},
  'p|pid=s'		=> \$opt{pid},
  'm|mapq=i'		=> \$opt{mapq},
  'q|mapqf=i'		=> \$opt{mapqf},
  'i|pidf=i'		=> \$opt{pidf},
  'g|gap=i'		=> \$opt{gap},
  'y|ci=i'		=> \$opt{ci},
  'sambamba=s'		=> \$opt{sambamba}
);

check_opt();

my %reads;
my %segments;
my %breakpoints;
my %breakpoints_region;
my %hanging_breakpoint_region;
my %structural_variants;
my %structural_variants_region;
my %structural_variants_region_2;

my $segmentID = 1;

my $readID = 1;
my $breakpointID = 1;
my $hanging_breakpointID = -1;
my $svID = 1;
my $sample = $1 if $opt{bam} =~ /(\w+)./;
my %svs;

my $date = strftime "%Y%m%d", localtime;

print_vcf_header();
parse_bam();
parse_reads();
parse_breakpoints();
parse_svs();

warn `date` . "Done\n";

sub parse_bam {
  warn `date` . "Busy with parsing bam file...\n";
  open BAM, "$opt{sambamba} view -t $opt{threads} $opt{bam} | ";
  while( <BAM> ) {
    chomp;    
    my @line = split/\t/, $_;
    my $read;
    if ( $reads{ $line[0] } ) {
      $read = $reads{ $line[0] };
    } else {
      $read = new Read( $line[0], $line[9], $line[5] );
      $reads{ $line[0] } = $read;
    }    
    next if $line[1] & 4 or $line[4] < $opt{mapq};
    my $segment = new Segment( $segmentID, $line[0], $line[1], $line[2], $line[3], $line[4], length( $line[9] ) );    
    $segment->parseCigar( $line[5] );
    next if $segment->{_pid} < $opt{pid};
    $read->addSegment( $segment );    
    $segments{ $segmentID } = $segment;
    $segmentID++;
  }
  close BAM;  
}

sub parse_reads { 
  warn `date` . "Busy with parsing read segments...\n";
  foreach my $qname ( keys %reads ) {
    my @clips = sort { $a <=> $b } keys %{ $reads{ $qname }->{_segments} };
    next if ( scalar( @clips ) == 1 or scalar( @clips ) > $opt{split} );
    foreach my $i ( 0..( $#clips - 1 ) ) {
      my $i2 = ( $i + 1 );      

      my $segment_1 = $segments{ $reads{ $qname }->{_segments}{ $clips[ $i ] } };
      my $segment_2 = $segments{ $reads{ $qname }->{_segments}{ $clips[ $i2 ] } };      
      $segment_1->setPlength( $reads{ $qname }->{_length} );
      $segment_2->setPlength( $reads{ $qname }->{_length} );
      
      my $breakpoint = new Breakpoint( $breakpointID++, $segment_1, $segment_2 );
                
      my $gap = ( $clips[ $i2 ] - ( $clips[ $i ] + $segment_1->{_length} ) );
      $breakpoint->setGap( $gap );     
      $breakpoint->setBreakpoint( $segment_1, $segment_2 );         
    
      if ( $segment_1->{_rname} gt $segment_2->{_rname} ) {
	$breakpoint->switchSegments();
      } elsif ( $segment_1->{_rname} =~ /^$segment_2->{_rname}$/ and $breakpoint->{_segment_1}{_pos} > $breakpoint->{_segment_2}{_pos}) {
	$breakpoint->switchSegments();
      }
    
      $breakpoint->setSVtype( );      
      $breakpoints{ $breakpoint->{_id} } = $breakpoint;
      $breakpoints_region{ join( "\t", $breakpoint->{_svtype}, $breakpoint->{_segment_1}{_rname}, $breakpoint->{_segment_2}{_rname}, $breakpoint->{_segment_1}{_flag}, $breakpoint->{_segment_2}{_flag} ) }{ $breakpoint->{_segment_1}{_pos} }{ $breakpoint->{_id} } = 1;
      
      if ( $i == 0 and $segment_1->{_clip} >= $opt{unmapped} ) {
	my $hanging_breakpoint_pos = $segment_1->{_pos};
	if ( $segment_1->{_flag} & 16 ) {
	  $hanging_breakpoint_pos = $segment_1->{_end};
	  $hanging_breakpoint_region{ $segment_1->{_rname} }{ 'T' }{ $hanging_breakpoint_pos }{ $hanging_breakpointID-- } = $segment_1->{_id};	  
	} else {
	  $hanging_breakpoint_region{ $segment_1->{_rname} }{ 'H' }{ $hanging_breakpoint_pos }{ $hanging_breakpointID-- } = $segment_1->{_id};
	}
      }
      if ( $i2 == $#clips and $segment_2->{_clip_2} >= $opt{unmapped} ) {
	my $hanging_breakpoint_pos = $segment_2->{_end};
	if ( $segment_2->{_flag} & 16 ) {	  
	  $hanging_breakpoint_pos = $segment_2->{_pos};
	  $hanging_breakpoint_region{ $segment_2->{_rname} }{ 'H' }{ $hanging_breakpoint_pos }{ $hanging_breakpointID-- } = $segment_2->{_id};
	} else {
	  $hanging_breakpoint_region{ $segment_2->{_rname} }{ 'T' }{ $hanging_breakpoint_pos }{ $hanging_breakpointID-- } = $segment_2->{_id};
	}
      }
    }
  }
}

sub parse_breakpoints {
  warn `date` . "Busy with parsing breakpoints...\n";
  foreach my $region ( keys %breakpoints_region ) {
    my $prev_pos_1 = -1;
    my %breakpoints_region_2;
    foreach my $pos_1 ( sort{ $a <=> $b } keys %{ $breakpoints_region{ $region } } ) {
      foreach my $breakpoint_id ( keys %{ $breakpoints_region{ $region }{ $pos_1 } } ) {
	if ( $prev_pos_1 == -1 ) {
	  $breakpoints_region_2{ $breakpoints{ $breakpoint_id }->{_segment_2}{_pos} }{ $pos_1 }{ $breakpoint_id } = 1;
	} elsif ( abs( $pos_1 - $prev_pos_1 ) <= $opt{distance} ) {
	  $breakpoints_region_2{ $breakpoints{ $breakpoint_id }->{_segment_2}{_pos} }{ $pos_1 }{ $breakpoint_id } = 1;
	} else {
	  parse_breakpoints_2( \%breakpoints_region_2 );
	  %breakpoints_region_2 = ();
	  $breakpoints_region_2{ $breakpoints{ $breakpoint_id }->{_segment_2}{_pos} }{ $pos_1 }{ $breakpoint_id } = 1;
	}
	$prev_pos_1 = $pos_1;
      }
    }
    parse_breakpoints_2( \%breakpoints_region_2 );
  }
}

sub parse_breakpoints_2 {
  my ( $breakpoints_region_2 ) = @_;
  my $prev_pos_1 = -1;
  my $prev_pos_2 = -1;
  my $sv;
  foreach my $pos_2 ( sort{ $a <=> $b } keys %{ $breakpoints_region_2 } ) {
    foreach my $pos_1 ( sort{ $a <=> $b } keys %{ $breakpoints_region_2->{ $pos_2 } } ) {
      foreach my $breakpoint_id ( keys %{ $breakpoints_region_2->{ $pos_2 }{ $pos_1 } } ) {
	my $breakpoint = $breakpoints{ $breakpoint_id };
	if ( $prev_pos_2 == -1 ) {
	  $sv = new SV( $svID++, $breakpoint );
	} elsif ( abs( $breakpoint->{_segment_2}{_pos} - $prev_pos_2 ) <= $opt{distance} ) {
	  $sv->addBreakpoint( $breakpoint );
# 	} elsif ( $breakpoint->{_segment_1}{_rname} =~ /^$breakpoint->{_segment_2}{_rname}$/ and $breakpoint->{_segment_1}{_pos} <= $prev_pos_2 and $breakpoint->{_segment_2}{_pos} >= $prev_pos_1 and abs( $breakpoint->{_segment_1}{_pos} - $prev_pos_1) <= $opt{distance} and abs( $breakpoint->{_segment_2}{_pos} - $prev_pos_2 ) <= $opt{distance} ) {
# 	} elsif ( $breakpoint->{_segment_1}{_rname} =~ /^$breakpoint->{_segment_2}{_rname}$/ and abs( $breakpoint->{_segment_2}{_pos} - $prev_pos_2 ) <= $opt{distance} ) {
# 	  $sv->addBreakpoint( $breakpoint );
# 	} elsif ( $breakpoint->{_segment_1}{_rname} !~ /^$breakpoint->{_segment_2}{_rname}$/ and abs( $breakpoint->{_segment_1}{_pos} - $prev_pos_1) <= $opt{distance} and abs( $breakpoint->{_segment_2}{_pos} - $prev_pos_2 ) <= $opt{distance} ) {
# 	} elsif ( $breakpoint->{_segment_1}{_rname} !~ /^$breakpoint->{_segment_2}{_rname}$/ and abs( $breakpoint->{_segment_2}{_pos} - $prev_pos_2 ) <= $opt{distance} ) {
# 	  $sv->addBreakpoint( $breakpoint );
	} else {    
	  if ( sum( @{ $sv->{_format}->{'DV'} } ) >= $opt{count}*2 ) {
	    if ( $sv->{_flag_1} & 16 ) {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_pos} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_pos} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[0]++;
		  $sv->{_format}{'VO'}->[0] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ [ $segments{ $segment_id }->{_pid} ], undef ] );
		  $sv->addInfoField( "MAPQ", [ [ $segments{ $segment_id }->{_mapq} ], undef ] );
		  $sv->addInfoField( "PLENGTH", [ [ $segments{ $segment_id }->{_plength} ], undef ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    } else {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_pos} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_pos} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[0]++;
		  $sv->{_format}{'VO'}->[0] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ [ $segments{ $segment_id }->{_pid} ], undef ] );
		  $sv->addInfoField( "MAPQ", [ [ $segments{ $segment_id }->{_mapq} ], undef ] );
		  $sv->addInfoField( "PLENGTH", [ [ $segments{ $segment_id }->{_plength} ], undef ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    }
	    if ( $sv->{_flag_2} & 16 ) {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'T' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_info}{'END'} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_info}{'END'} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'T' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[1]++;
		  $sv->{_format}{'VO'}->[1] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ undef, [ $segments{ $segment_id }->{_pid} ] ] );
		  $sv->addInfoField( "MAPQ", [ undef, [ $segments{ $segment_id }->{_mapq} ] ] );
		  $sv->addInfoField( "PLENGTH", [ undef, [ $segments{ $segment_id }->{_plength} ] ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    } else {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'H' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_pos} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_pos} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'H' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[1]++;
		  $sv->{_format}{'VO'}->[1] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ undef, [ $segments{ $segment_id }->{_pid} ] ] );
		  $sv->addInfoField( "MAPQ", [ undef, [ $segments{ $segment_id }->{_mapq} ] ] );
		  $sv->addInfoField( "PLENGTH", [ undef, [ $segments{ $segment_id }->{_plength} ] ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    }
	    $structural_variants{ $sv->{_id} } = $sv;
	    $structural_variants_region{ $sv->{_chr} }{ min( @{ $sv->{_pos} } ) }{ max( @{ $sv->{_pos} } ) }{ $sv->{_id} } = 0;
	    $structural_variants_region{ $sv->{_chr2} }{ min( @{ $sv->{_info}{'END'} } ) }{ max( @{ $sv->{_info}{'END'} } ) }{ $sv->{_id} } = 1;
	    $structural_variants_region_2{ $sv->{_chr} }{ $sv->{_id } } = 1;

	  }
	  $sv = new SV( $svID++, $breakpoint );
	}
	$sv->addInfoField( "PID", [ [ $segments{ $breakpoint->{_segment_1}{_id} }->{_pid} ], [ $segments{ $breakpoint->{_segment_2}{_id} }->{_pid} ] ] );
	$sv->addInfoField( "MAPQ", [ [ $segments{ $breakpoint->{_segment_1}{_id} }->{_mapq} ], [ $segments{ $breakpoint->{_segment_2}{_id} }->{_mapq} ] ] );
	$sv->addInfoField( "PLENGTH", [ [ $segments{ $breakpoint->{_segment_1}{_id} }->{_plength} ], [ $segments{ $breakpoint->{_segment_2}{_id} }->{_plength} ] ] );
	
	$sv->addInfoField( "RLENGTH", [ $reads{ $segments{$breakpoint->{_segment_2}{_id}}->{_qname} }->{_length} ] );
	$sv->addInfoField( "GAP", [ $breakpoint->{_gap} ] );
	if ( $segments{$breakpoint->{_segment_2}{_id}}->{_qname} =~ /2D_2d$/ ) {
	  $sv->addInfoField( "RT", [ "2d" ] );
	} elsif ( $segments{$breakpoint->{_segment_2}{_id}}->{_qname} =~ /2D_complement$/ ) {
	  $sv->addInfoField( "RT", [ "complement" ] );
	} else {
	  $sv->addInfoField( "RT", [ "template" ] );
	}

	$prev_pos_1 = $breakpoint->{_segment_1}{_pos};
	$prev_pos_2 = $pos_2;
      }
    }
  }
  if ( sum( @{ $sv->{_format}->{'DV'} } ) >= $opt{count}*2 ) {
    if ( $sv->{_flag_1} & 16 ) {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_pos} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_pos} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[0]++;
		  $sv->{_format}{'VO'}->[0] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ [ $segments{ $segment_id }->{_pid} ], undef ] );
		  $sv->addInfoField( "MAPQ", [ [ $segments{ $segment_id }->{_mapq} ], undef ] );
		  $sv->addInfoField( "PLENGTH", [ [ $segments{ $segment_id }->{_plength} ], undef ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    } else {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_pos} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_pos} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[0]++;
		  $sv->{_format}{'VO'}->[0] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ [ $segments{ $segment_id }->{_pid} ], undef ] );
		  $sv->addInfoField( "MAPQ", [ [ $segments{ $segment_id }->{_mapq} ], undef ] );
		  $sv->addInfoField( "PLENGTH", [ [ $segments{ $segment_id }->{_plength} ], undef ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    }
	    if ( $sv->{_flag_2} & 16 ) {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'T' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_info}{'END'} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_info}{'END'} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'T' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'T' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[1]++;
		  $sv->{_format}{'VO'}->[1] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ undef, [ $segments{ $segment_id }->{_pid} ] ] );
		  $sv->addInfoField( "MAPQ", [ undef, [ $segments{ $segment_id }->{_mapq} ] ] );
		  $sv->addInfoField( "PLENGTH", [ undef, [ $segments{ $segment_id }->{_plength} ] ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    } else {
	      foreach my $hanging_pos ( sort{ $a <=> $b } keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'H' } } ) {
		next if $hanging_pos < ( min( @{ $sv->{_pos} } ) - $opt{distance} );
		last if $hanging_pos > ( max( @{ $sv->{_pos} } ) + $opt{distance} );
		foreach my $hanging_id ( keys %{ $hanging_breakpoint_region{ $sv->{_chr2} }{ 'H' }{ $hanging_pos } } ) {
		  my $segment_id = $hanging_breakpoint_region{ $sv->{_chr} }{ 'H' }{ $hanging_pos }{ $hanging_id };
		  $sv->{_format}{'HR'}->[1]++;
		  $sv->{_format}{'VO'}->[1] += ( 1-10**( -$segments{ $segment_id }->{_mapq}/10.0 ) );
		  $sv->addInfoField( "PID", [ undef, [ $segments{ $segment_id }->{_pid} ] ] );
		  $sv->addInfoField( "MAPQ", [ undef, [ $segments{ $segment_id }->{_mapq} ] ] );
		  $sv->addInfoField( "PLENGTH", [ undef, [ $segments{ $segment_id }->{_plength} ] ] );
		  $sv->addInfoField( "RLENGTH", [ $reads{ $segments{ $segment_id }->{_qname} }->{_length} ] );
		  if ( $segments{ $segment_id }->{_qname} =~ /2D_2d$/ ) {
		    $sv->addInfoField( "RT", [ "2d" ] );
		  } elsif ( $segments{ $segment_id }->{_qname} =~ /2D_complement$/ ) {
		    $sv->addInfoField( "RT", [ "complement" ] );
		  } else {
		    $sv->addInfoField( "RT", [ "template" ] );
		  }
		}
	      }
	    }
    $structural_variants{ $sv->{_id} } = $sv;
    $structural_variants_region{ $sv->{_chr} }{ min( @{ $sv->{_pos} } ) }{ max( @{ $sv->{_pos} } ) }{ $sv->{_id} } = 0;
    $structural_variants_region{ $sv->{_chr2} }{ min( @{ $sv->{_info}{'END'} } ) }{ max( @{ $sv->{_info}{'END'} } ) }{ $sv->{_id} } = 1;
    $structural_variants_region_2{ $sv->{_chr} }{ $sv->{_id } } = 1;
  }
}

sub parse_svs {
  warn `date` .  "Busy with reference reads...\n";
  my $prev_rname = 0;
  foreach my $segment_id ( sort{ $a <=> $b } keys %segments ) {
    if ( $segments{ $segment_id }->{_rname} ne $prev_rname and $prev_rname != 0 ) {
      delete $structural_variants_region{ $prev_rname };
    }
    foreach my $sv_min ( sort{ $a <=> $b} keys %{ $structural_variants_region{ $segments{ $segment_id }->{_rname} } } ) {
      if ( $sv_min <= $segments{ $segment_id }->{_pos} ) {
	delete $structural_variants_region{ $segments{ $segment_id }->{_rname} }{ $sv_min };
	next;
      }
      last if $sv_min > $segments{ $segment_id }->{_end};
      foreach my $sv_max ( sort{ $a <=> $b} keys %{ $structural_variants_region{ $segments{ $segment_id }->{_rname} }{ $sv_min } } ) {
	last if $sv_max >= $segments{ $segment_id }->{_end};
	foreach my $sv_id ( keys %{ $structural_variants_region{ $segments{ $segment_id }->{_rname} }{ $sv_min }{ $sv_max } } ) {
	  my $x = $structural_variants_region{ $segments{ $segment_id }->{_rname} }{ $sv_min }{ $sv_max }{ $sv_id };
	  if ( $sv_min > ( $segments{ $segment_id }{_pos} + $opt{refdistance} ) and $sv_max < ( $segments{ $segment_id }{_end} - $opt{refdistance} ) ) {
	    $structural_variants{ $sv_id }->{_format}{'DR'}->[$x]++;
	    $structural_variants{ $sv_id }->{_format}{'RO'}->[$x] += ( 1-10**( -$segments{ $segment_id }{_mapq}/10.0 ) );
	  }
	}
      }
    }
    $prev_rname = $segments{ $segment_id }->{_rname};
  }
  warn `date` .  "Busy with printing to vcf...\n";
  foreach my $sv_id ( sort{ $a <=> $b } keys %structural_variants ) {    
    $structural_variants{ $sv_id }->setArguments() unless $structural_variants{ $sv_id }->{_set} == 1;   
    foreach my $sv_id_2 ( sort{ $a <=> $b } keys %{ $structural_variants_region_2{ $structural_variants{ $sv_id }->{_chr} } } ) {
      next unless $sv_id_2 > $sv_id;
      $structural_variants{ $sv_id_2 }->setArguments() unless $structural_variants{ $sv_id_2 }->{_set} == 1;      
      next unless $structural_variants{ $sv_id }->{_chr} eq $structural_variants{ $sv_id_2 }->{_chr};
      next unless $structural_variants{ $sv_id }->{_chr2} eq $structural_variants{ $sv_id_2 }->{_chr2};
      if ( abs( $structural_variants{ $sv_id }->{_pos} - $structural_variants{ $sv_id_2 }->{_pos} ) <= $opt{window} and abs( $structural_variants{ $sv_id }->{_info}{'END'} - $structural_variants{ $sv_id_2 }->{_info}{'END'} ) <= $opt{window} ) {
	
	$structural_variants{ $sv_id }{'cluster'}++;
	$structural_variants{ $sv_id_2 }{'cluster'}++;
      }
      next unless $structural_variants{ $sv_id_2 }->{_info}{'SVTYPE'} eq "BND";
      next unless $structural_variants{ $sv_id }->{_info}{'SVTYPE'} eq "BND";
      next if $structural_variants{ $sv_id }->{_alt} =~ /\](\w+):(\d+)\]\w+/ and $structural_variants{ $sv_id_2 }->{_alt} !~ /\w+\[(\w+):(\d+)\[/;
      next if $structural_variants{ $sv_id }->{_alt} =~ /\[(\w+):(\d+)\[\w+/ and $structural_variants{ $sv_id_2 }->{_alt} !~ /\w+\](\w+):(\d+)\]/;
      next if $structural_variants{ $sv_id }->{_alt} =~ /\w+\](\w+):(\d+)\]/ and $structural_variants{ $sv_id_2 }->{_alt} !~ /\[(\w+):(\d+)\[\w+/;
      next if $structural_variants{ $sv_id }->{_alt} =~ /\w+\[(\w+):(\d+)\[/ and $structural_variants{ $sv_id_2 }->{_alt} !~ /\](\w+):(\d+)\]\w+/;
      next unless abs( $structural_variants{ $sv_id }->{_pos} - $structural_variants{ $sv_id_2 }->{_pos} ) <= $opt{matedistance};
      next unless abs( $structural_variants{ $sv_id }->{_info}{'END'} - $structural_variants{ $sv_id_2 }->{_info}{'END'} ) <= $opt{matedistance};
      $structural_variants{ $sv_id }->{_info}{'MATEID'} = $sv_id_2;
      $structural_variants{ $sv_id_2 }->{_info}{'MATEID'} = $sv_id;
    }
    push @{ $structural_variants{ $sv_id }->{_filter} }, "SVcluster" if $structural_variants{ $sv_id }{'cluster'} > $opt{cluster};
    push @{ $structural_variants{ $sv_id }->{_filter} }, "GAP" if $structural_variants{ $sv_id }->{_info}{'GAP'} > $opt{gap} and $structural_variants{ $sv_id }->{_info}{'SVTYPE'} ne "INS";    
    if ( $structural_variants{ $sv_id }->{_info}{'MAPQ'} =~ /(\d+),(\d+)/ ) {
      push @{ $structural_variants{ $sv_id }->{_filter} }, "MapQual" if $1 < $opt{mapqf} or $2 < $opt{mapqf};
    }
    if ( $structural_variants{ $sv_id }->{_info}{'PID'} =~ /(\d.\d+),(\d.\d+)/ ) {
      push @{ $structural_variants{ $sv_id }->{_filter} }, "PID" if $1 < $opt{pidf} or $2 < $opt{pidf};
    }
    if ( $structural_variants{ $sv_id }->{_info}{'CIPOS'} =~ /(\d+),(\d+)/ ) {
      push @{ $structural_variants{ $sv_id }->{_filter} }, "CIPOS" if $1 > $opt{ci} or $2 > $opt{ci};
    }
    if ( $structural_variants{ $sv_id }->{_info}{'CIEND'} =~ /(\d+),(\d+)/ ) {
      push @{ $structural_variants{ $sv_id }->{_filter} }, "CIEND" if $1 > $opt{ci} or $2 > $opt{ci};
    }
    $structural_variants{ $sv_id }->printVCF();
  }
}

sub check_opt {  
  die usage() if $opt{help};
  die usage("No bam file given") unless $ARGV[$#ARGV] =~ /.bam$/;
  die usage("Bam file not found") unless -e $ARGV[$#ARGV];
  my $sambamba_check = `$opt{sambamba} 2>&1 1>/dev/null`;
  die usage("sambamba not found") unless $sambamba_check =~ /Usage: sambamba \[command\]/;
  $opt{bam} = $ARGV[$#ARGV];
}

sub print_vcf_header {
  print <<END;
##fileformat=VCFv4.1
##fileDate=$date
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=BND,Description="Breakend">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Distance between the two genomic positions ( END - POS )">
##INFO=<ID=RT,Number=3,Type=Integer,Description="Number of the different read types ( 2d , template , complement )">
##INFO=<ID=GAP,Number=1,Type=Integer,Description="Median number of bases between the two segments of the SV, in case of an insertion this is the size of the insertion">
##INFO=<ID=MAPQ,Number=2,Type=Integer,Description="Median mapping quality of the two segments of the structural variant">
##INFO=<ID=PID,Number=2,Type=Float,Description="Median percentage identity to the reference of the two segments of the structural variant">
##INFO=<ID=PLENGTH,Number=2,Type=Float,Description="Median segment length percentage of the two segments of the structural variant">
##INFO=<ID=RLENGTH,Number=1,Type=Integer,Description="Median length of the total reads">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Phred-scaled probability that this site is variant (non-reference in this sample)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HR,Number=2,Type=Integer,Description="Number of hanging variant reads">
##FILTER=<ID=SVcluster,Description="There are more than $opt{cluster} SVs in a window of $opt{window} on both sides">
##FILTER=<ID=GAP,Description="The median gap size is larger than $opt{gap} for non insertions">
##FILTER=<ID=MapQual,Description="The median mapping quality is less than $opt{mapqf}">
##FILTER=<ID=PID,Description="The PID of one of the segments is less than $opt{pidf}">
##FILTER=<ID=CIPOS,Description="The CIPOS is greater or less than $opt{ci}">
##FILTER=<ID=CIEND,Description="The CIEND is greater or less than $opt{ci}">
END
  print join( "\t", qw( #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ), $sample ) . "\n";  
}
