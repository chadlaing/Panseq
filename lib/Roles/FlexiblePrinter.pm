#!/usr/bin/env perl

package Roles::FlexiblePrinter;
#allows a unified, flexible printing strategy
use strict;
use warnings;
use FindBin;
use IO::File;
use lib "$FindBin::Bin/../";
use Role::Tiny;

sub outputFH{
	my $self=shift;
	$self->{'_outputFH'}=shift // return $self->{'_outputFH'};
}

#methods
sub printOut{
	my($self)=shift;
	
	my $outFH = $self->outputFH // *STDOUT{IO};
	print $outFH @_;
}

sub DESTROY{
	my $self=shift;

	if($self->outputFH){
		$self->outputFH->close();
	}
}

1;
