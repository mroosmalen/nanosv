#! /usr/bin/perl

package Breakpoint;

use strict;

sub new {
  my ( $class, $id, $segment_1, $segment_2 ) = @_;
  my %self = (
    '_id' => $id,
    '_segment_1' => { 
	'_id' => $segment_1->{_id},
	'_rname' => $segment_1->{_rname},
	'_flag' => $segment_1->{_flag},
	'_mapq' => $segment_1->{_mapq}
    },
    '_segment_2' => {
	'_id' => $segment_2->{_id},
	'_rname' => $segment_2->{_rname},
	'_flag' => $segment_2->{_flag},
	'_mapq' => $segment_2->{_mapq}
    }
  );
  
  return bless \%self, $class;
}

sub setGap {
  my ( $self, $gap ) = @_;
  $self->{_gap} = $gap;
}

sub setBreakpoint {
  my ( $self, $segment_1, $segment_2 ) = @_;
  $self->{_segment_1}{_pos} = $segment_1->{_flag} & 16 ? $segment_1->{_pos} : $segment_1->{_end};
  $self->{_segment_2}{_pos} = $segment_2->{_flag} & 16 ? $segment_2->{_end} : $segment_2->{_pos};
}

sub setPlength {
  my ( $self, $segment_1, $segment_2, $read ) = @_;
  $self->{_segment_1}{_plength} = ( $segment_1->{_length} / $read->{_length} );
  $self->{_segment_2}{_plength} = ( $segment_2->{_length} / $read->{_length} );
}

sub switchSegments {
  my ( $self ) = @_;
  ( $self->{_segment_1}, $self->{_segment_2} ) = ( $self->{_segment_2}, $self->{_segment_1} );
  $self->{_segment_1}{_flag} = $self->{_segment_1}{_flag} & 16 ? 0 : 16;
  $self->{_segment_2}{_flag} = $self->{_segment_2}{_flag} & 16 ? 0 : 16;  
}

sub setSVtype {
  my ( $self ) = @_;
  if ( $self->{_segment_1}{_rname} !~ /^$self->{_segment_2}{_rname}$/ ) {
      $self->{_svtype} = "BND";
  } else {
    $self->{_svtype} = $self->{_segment_1}{_flag} & 16 ? $self->{_segment_2}{_flag} & 16 ? "DUP" : "BND" : $self->{_segment_2}{_flag} & 16 ? "BND" : "DEL";
#     if ( abs( $self->{_segment_2}{_pos} - $self->{_segment_1}{_pos} ) < $self->{_gap} and $self->{_svtype} eq "DEL") {
#       $self->{_svtype} = "INS";
#     }
    if ( abs( $self->{_segment_2}{_pos} - $self->{_segment_1}{_pos} ) < $self->{_gap} ) {
      $self->{_svtype} = "INS";
      $self->{_segment_2}{_pos} = $self->{_segment_1}{_pos}+1;
      $self->{_segment_1}{_flag} = 0;
      $self->{_segment_2}{_flag} = 0;
    }
    
  }
}

1;