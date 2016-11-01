#! /usr/bin/perl

package SV;

use strict;
use List::Util qw( min max sum );
use Data::Dumper;

sub new {
  my ( $class, $id, $breakpoint ) = @_;
  
  my %self = (
    '_id' => $id,
    '_chr' => $breakpoint->{_segment_1}{_rname},
    '_chr2' => $breakpoint->{_segment_2}{_rname},
    '_flag_1' => $breakpoint->{_segment_1}{_flag},
    '_flag_2' => $breakpoint->{_segment_2}{_flag},
    '_pos' => [ $breakpoint->{_segment_1}{_pos} ],
    '_ref' => 'N',
    '_alt' => "<$breakpoint->{_svtype}>",
    '_qual' => '.',
    '_filter' => [ 'PASS' ],
    '_info' => {
      'PRECISE' => "IMPRECISE",
      'END' => [ $breakpoint->{_segment_2}{_pos} ],
      'SVTYPE' => $breakpoint->{_svtype},
      'SVMETHOD' => 'nanosv',
    },
    '_format' => {
      'GT' => './.',
      'DV' => [ 1, 1],
      'VO' => [ ( 1-10**( -$breakpoint->{_segment_1}{_mapq}/10.0 ) ), ( 1-10**( -$breakpoint->{_segment_2}{_mapq}/10.0 ) ) ],
      'DR' => [ 0, 0 ],
      'RO' => [ 0, 0 ],
      'HR' => [ 0, 0 ]
    },      
    '_breakpoints' => [ $breakpoint->{_id} ],
    '_set' => 0
  );
  
  return bless \%self, $class;
}

sub addBreakpoint {
  my ( $self, $breakpoint ) = @_;
  push @{ $self->{_breakpoints} }, $breakpoint->{_id};
  push @{ $self->{_pos} }, $breakpoint->{_segment_1}{_pos};
  push @{ $self->{_info}{'END'} }, $breakpoint->{_segment_2}{_pos};
  $self->{_format}{'DV'}->[0]++;
  $self->{_format}{'DV'}->[1]++;
  $self->{_format}{'VO'}->[0] += ( 1-10**( -$breakpoint->{_segment_1}{_mapq}/10.0 ) );
  $self->{_format}{'VO'}->[1] += ( 1-10**( -$breakpoint->{_segment_2}{_mapq}/10.0 ) );
}

sub addInfoField {
  my ( $self, $key, $value ) = @_;
  if ( $self->{_info}{$key} =~ /ARRAY/ ) {
    if ( $self->{_info}{$key}->[0] =~ /ARRAY/ or $self->{_info}{$key}->[1] =~ /ARRAY/ ) {
      push @{ $self->{_info}{$key}->[0] }, $value->[0][0] if $value->[0];
      push @{ $self->{_info}{$key}->[1] }, $value->[1][0] if $value->[1];
    } else {
      push @{ $self->{_info}{$key} }, $value->[0];
    }
  } else {
    $self->{_info}{$key} = $value;  
  }
}

sub addFormatField {
  my ( $self, $key, $value ) = @_;
  $self->{_format}{$key} = $value;
}

sub setArguments {
  my ( $self ) = @_;
 
#   warn Dumper($self) . "\n";
 
  $self->{_info}{'CIPOS'} = ( min( @{ $self->{_pos} } ) - median( @{ $self->{_pos} } ) ) . "," . ( max( @{ $self->{_pos} } ) - median( @{ $self->{_pos} } ) );
  $self->{_info}{'CIEND'} = ( min( @{ $self->{_info}{'END'} } ) - median( @{ $self->{_info}{'END'} } ) ) . "," . ( max( @{ $self->{_info}{'END'} } ) - median( @{ $self->{_info}{'END'} } ) );
  $self->{_info}{'PRECISE'} = "PRECISE" if $self->{_info}{'CIPOS'} eq "0,0" and $self->{_info}{'CIEND'};
  $self->{_pos} = median( @{ $self->{_pos} } );
  $self->{_info}{'END'} = median( @{ $self->{_info}{'END'} } );  
  $self->{_info}{'SVLEN'} = ( $self->{_info}{'END'} - $self->{_pos} );    
  $self->setInfoField();
  
  my $dup = 0;
  $dup = 1  if $self->{_info}{'SVTYPE'} eq "DUP";
  
  my @gt_lplist = bayes_gt( sum( @{ $self->{_format}{'RO'} } ), sum( @{ $self->{_format}{'VO'} } ), $dup );
  my $gt_idx = 0;
  $gt_lplist[$gt_idx] > $gt_lplist[$_] or $gt_idx = $_ for 1 .. $#gt_lplist;
     
  my $gt_sum = 0;
  foreach my $gt (@gt_lplist) {
    $gt_sum += 10**$gt;
  }
  
  if ( $gt_sum > 0 ) {
    my $gt_sum_log = log10( $gt_sum );
    my $sample_qual = abs(-10 * ( $gt_lplist[0] - $gt_sum_log ) );
    my $phred_gq = -1;
    if ( 1 - ( 10**$gt_lplist[$gt_idx] / 10**$gt_sum_log ) == 0 ) {
      $phred_gq = 200;
    } else {
      $phred_gq = abs(-10 * log10( 1 - ( 10**$gt_lplist[ $gt_idx ] / 10**$gt_sum_log ) ) );
    }
    $self->{_format}{'GQ'} = int( $phred_gq );
    $self->{_format}{'SQ'} = sprintf( "%.3f", $sample_qual );
    
    if ( $gt_idx == 0 ) {
      $self->{_format}{'GT'} = '0/0';
    } elsif ( $gt_idx == 1 ) {
      $self->{_format}{'GT'} = '0/1';
    } elsif ( $gt_idx == 2 ) {
      $self->{_format}{'GT'} = '1/1';
    } 
  }
  
  if ( $self->{_alt} eq "<BND>" ) {
    if ( $self->{_flag_1} & 16 ) {
      if ( $self->{_flag_2} & 16 ) {
	$self->{_alt} = "]" . $self->{_chr2} . ":" . $self->{_info}{'END'} . "]" . $self->{_ref};
      } else {
	$self->{_alt} = "[" . $self->{_chr2} . ":" . $self->{_info}{'END'} . "[" . $self->{_ref};
      }
    } else {
      if ( $self->{_flag_2} & 16 ) {
	$self->{_alt} = $self->{_ref} . "]" . $self->{_chr2} . ":" . $self->{_info}{'END'} . "]";
      } else {
	$self->{_alt} = $self->{_ref} . "[" . $self->{_chr2} . ":" . $self->{_info}{'END'} . "[";
      }
    }
  }
  
  $self->{_set} = 1;
  
}

sub bayes_gt {
  my ( $ref, $alt, $dup ) = @_;
  my @p_alt = ( 0.01, 0.5, 0.9 );
  @p_alt = ( 0.01, 0.3, 0.5 ) if $dup;
  my $total = ( $ref + $alt );

  my $lp_homref = log_choose( $total, $alt ) + $alt * log10( $p_alt[0] ) + $ref * log10( 1 - $p_alt[0] );
  my $lp_het = log_choose( $total, $alt ) + $alt * log10( $p_alt[1] ) + $ref * log10( 1 - $p_alt[1] );
  my $lp_homalt = log_choose( $total, $alt) + $alt * log10( $p_alt[2] ) + $ref * log10( 1 - $p_alt[2] );

  return( $lp_homref, $lp_het, $lp_homalt );
}

sub log_choose {
  my ( $n, $k ) = @_;
  my $r = 0.0;
  if ( $k * 2 > $n ) {
    $k = $n - $k;
  }
  
  foreach my $d ( 1..$k ) {
    $r += log10($n);
    $r -= log10($d);
    $n -= 1;
  }
  
  return $r;
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub printVCF {
  my ( $self ) = @_;  

  if ( scalar @{ $self->{_filter} } == 1 ) {
    $self->{_filter} = "PASS";
  } else {
    $self->{_filter} = join(",", @{ $self->{_filter} } );
    $self->{_filter} =~ s/PASS,//;
  }
  
  print $self->{_chr} . "\t";
  print $self->{_pos} . "\t";
  print $self->{_id} . "\t";
  print $self->{_ref} . "\t";
  print $self->{_alt} . "\t";
  print $self->{_qual} . "\t";
  print $self->{_filter} . "\t";
  print $self->{_info}{'PRECISE'};
  foreach my $field ( keys %{ $self->{_info} } ) {
    next if $field eq 'PRECISE';
    next if $field eq "END" and $self->{_info}{'SVTYPE'} eq "BND";
    next if $field eq "SVLEN" and $self->{_chr} !~ /^$self->{_chr2}$/;
    print ";" . $field . "=" . $self->{_info}{ $field };
  }
  print "\t" . 'GT';
  my @values;
  foreach my $field ( keys %{ $self->{_format} } ) {
    next if $field eq 'GT';
    next unless $field =~ /DR|DV|HR|SQ|GQ/;
#     next if $field eq "SQ" or $field eq "GQ" or $field eq "SR" or $field eq "SV";
    print ":" . $field;
    my $value = $self->{_format}{ $field };
    $value = join(",",@{ $self->{_format}{ $field } } ) if $value =~ /ARRAY/;
    push @values, $value;
  }
  print "\t" . $self->{_format}{'GT'} . ":" . join( ":", @values ) . "\n";
}

sub setInfoField {
  my ( $self ) = @_;
  foreach my $field (keys %{$self->{_info}} ) {
    if ( $field eq "RT" ) {
      my @rt = ( 0, 0, 0 );
      foreach my $type ( @{ $self->{_info}{ $field } } ) {
	$rt[0]++ if $type eq "2d";
	$rt[1]++ if $type eq "template";
	$rt[2]++ if $type eq "complement";
      }
      $self->{_info}{ $field } = join( ",", @rt );
    } elsif ( $self->{_info}{ $field } =~ /ARRAY/ ) {
      if ( $self->{_info}{ $field }->[0] =~ /ARRAY/ ) {
	$self->{_info}{ $field }->[0] = median( @{ $self->{_info}{ $field }->[0] } );
	$self->{_info}{ $field }->[1] = median( @{ $self->{_info}{ $field }->[1] } );
	$self->{_info}{ $field } = join(",", @{ $self->{_info}{ $field } } );
      } else {
	$self->{_info}{ $field } = median( @{ $self->{_info}{ $field } } );
      }       
    }
  }
}

sub median {
  my ( @array ) = @_;
  if ( $array[0] =~ /\./ ) {
    return sprintf( "%.3f", sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), int( $#_/2 + 0.5 ) ] )/2 );
  } else {
    return int ( sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), int( $#_/2 + 0.5 ) ] )/2 );
  }
}

1;