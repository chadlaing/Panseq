#!/usr/bin/perl

package DeltaBlockObject;
use Carp;

use Object::Tiny::RW qw{
  refName
  queryName
  refLength
  queryLength
  refStart
  refEnd
  queryStart
  queryEnd
  gapArray
  errors
  similarityErrors
  stopCodons
};

#methods
#given a block returns its approximate percent ID corresponding to the length of query sequence
sub getPercentID {
	my $self = shift;

	my $gapsQuery          = $self->getGapNumQuery();
	my $refAlignmentLength = $self->refEnd() - $self->refStart() + 1;
	my $totalErrors =
	  $self->errors() +
	  ( $self->queryLength() - $refAlignmentLength ) +
	  $gapsQuery;
	  
	my $percID =
	  ( $self->queryLength() + $gapsQuery - $totalErrors ) /
	  ( $self->queryLength() + $gapsQuery );

	return $percID;
}

sub getGapNumQuery {
	my $self      = shift;
	my $gapsQuery = 0;
	foreach ( @{ $self->gapArray() } ) {
		if ( $_ > 0 ) {
			$gapsQuery++;
		}
	}
	my $queryAlignmentLength;
	if ( $self->queryEnd() > $self->queryStart() ) {
		$queryAlignmentLength = $self->queryEnd() - $self->queryStart() + 1;
	}
	else {
		$queryAlignmentLength =
		  $self->queryStart() - $self->queryEnd() + 1;
	}
	$gapsQuery += ( $self->queryLength() - $queryAlignmentLength );

	return $gapsQuery;
}

sub setAbsoluteStartEndValues{
	my($self)=shift;
	
	if($self->refStart > $self->refEnd){
		my $temp = $self->refStart;
		$self->refStart($self->refEnd);
		$self->refEnd($temp);
	}
	
	if($self->queryStart > $self->queryEnd){
		my $temp = $self->queryStart;
		$self->queryStart($self->queryEnd);
		$self->queryEnd($temp);
	}
}
1;
