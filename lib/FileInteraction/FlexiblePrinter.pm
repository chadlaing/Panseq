#!/usr/bin/perl

package FileInteraction::FlexiblePrinter;
#allows a unified, flexible printing strategy
use strict;
use warnings;
use diagnostics;
use FindBin;
use lib "$FindBin::Bin";

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_initialize(@_);
    return $self;
}


#get/set
sub _initialize{
	my $self=shift;
	
	my %init = @_;	
	$self->outputFH($init{'outputFH'}) if defined $init{'outputFH'};
}
 
sub outputFH{
	my($self)=shift;
	
	if(@_){
		$self->{'_outputFilehandle'}=shift;
	}
	else{
		return $self->{'_outputFilehandle'} || *STDOUT{IO};
	}
}

#methods
sub print{
	my($self)=shift;
	
	my $outFile = $self->outputFH;
	print $outFile @_;
}

1;
