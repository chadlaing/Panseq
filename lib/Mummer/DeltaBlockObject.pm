#!/usr/bin/perl

package Mummer::DeltaBlockObject;
use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin";

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

sub stopCodons{
	my $self=shift;
	$self->{'_stopCodons'}=shift // return $self->{'_stopCodons'};
}

sub similarityErrors{
	my $self=shift;
	$self->{'_similarityErrors'}=shift // return $self->{'_similarityErrors'};
}

sub errors{
	my $self=shift;
	$self->{'_errors'}=shift // return $self->{'_errors'};
}

sub gapArray{
	my $self=shift;
	$self->{'_gapArray'}=shift // return $self->{'_gapArray'};
}

sub queryEnd{
	my $self=shift;
	$self->{'_queryEnd'}=shift // return $self->{'_queryEnd'};
}

sub queryStart{
	my $self=shift;
	$self->{'_queryStart'}=shift // return $self->{'_queryStart'};
}

sub refEnd{
	my $self=shift;
	$self->{'_refEnd'}=shift // return $self->{'_refEnd'};
}

sub refStart{
	my $self=shift;
	$self->{'_refStart'}=shift // return $self->{'_refStart'};
}

sub queryLength{
	my $self=shift;
	$self->{'_queryLength'}=shift // return $self->{'_queryLength'};
}

sub refLength{
	my $self=shift;
	$self->{'_refLength'}=shift // return $self->{'_refLength'};
}

sub queryName{
	my $self=shift;
	$self->{'_queryName'}=shift // return $self->{'_queryName'};
}

sub refName{
	my $self=shift;
	$self->{'_refName'}=shift // return $self->{'_refName'};
}


#methods
#given a block returns its approximate percent ID corresponding to the length of query sequence

sub _initialize{
	my $self=shift;
}


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
