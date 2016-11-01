#! /usr/bin/perl

package Segment;

use strict;
use Data::Dumper;

sub new {
  my ( $class, $id, $qname, $flag, $rname, $pos, $mapq, $length ) = @_;  
  my %self = (
    '_id' => $id,
    '_qname' => $qname,
    '_flag' => $flag,
    '_rname' => $rname,
    '_pos' => $pos,
    '_mapq' => $mapq,
    '_end' => ( $pos + $length ),
    '_length' => $length,
    '_clip' => undef,
    '_pid' => undef
  );
  
  return bless \%self, $class;
}

sub parseCigar {
  my ( $self, $cigar ) = @_;
  
  if ( $self->{_flag} & 16 ) {
    if ( $cigar =~ /(\d+)[HS]$/ ) {
      $self->{_clip} = $1;
    }
  } else {
    if ( $cigar =~ /^(\d+)[HS]/ ) {
      $self->{_clip} = $1;
    }
  }  
  
  my %regex;
  $regex{$2} += $1 while $cigar =~ /(\d+)([DI=])/g;  
  $self->{_end} += $regex{'D'};
  $self->{_end} -= $regex{'I'};
  $self->{_pid} = sprintf( "%.3f", ( $regex{'='} / $self->{_length} ) );
}

sub setBreakpoint {
  my ( $self, $index ) = @_;
  if ( $index == 1 ) {
    $self->{_breakpoint} = $self->{_flag} & 16 ? $self->{_pos} : $self->{_end};
  } elsif ( $index == 2 ) {
    $self->{_breakpoint} = $self->{_flag} & 16 ? $self->{_end} : $self->{_pos};
  }
}

sub setPlength {
  my ( $self, $rlength ) = @_;
  $self->{_plength} = sprintf( "%.3f", ( $self->{_length} / $rlength ) );
}

1;