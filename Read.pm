#! /usr/bin/perl

package Read;

use strict;

sub new {
  my ( $class, $qname, $seq, $cigar ) = @_;  

  my $left_length = $1 if $cigar =~ /^(\d+)H/;
  my $right_length = $1 if $cigar =~ /(\d+)H$/;
  
  my %self = (
    '_qname' => $qname,
    '_length' => ( length( $seq ) + $left_length + $right_length ),
    '_segments' => { }
  );  
  return bless \%self, $class;
}

sub addReadLength {
  my ( $self, $length ) = @_;
  $self->{_length} += $length;
}

sub addSegment {
  my ( $self, $segment ) = @_;
  $self->{_segments}->{ $segment->{_clip} } = $segment->{_id};
}

1;

