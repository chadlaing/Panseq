#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";

=head1 Description

Takes the output tables of Panseq (core_snps.txt) and (pan_genome.txt) and produces new tables
based on the specified parameters.
For example, take only SNP lines that contain fewer than X '-', and/or that have
at least 2 SNP characters present at least Y times.
Or, for pan_genome data, produce a binary table based on a new %ID cutoff value.

=cut

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


=head3 _initialize

Initializes the logger.
Assigns all values to class variables.
Anything else that the _initialize function does.

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 
    $self->logger->debug("Logger initialized in Modules::PanGenome::TableExtractor");  

    my %params = @_;

    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::PanGenome::TableExtractor");
		}
	}	
}

sub inputFile{
	my $self=shift;
	$self->{'_inputFile'}=shift // return $self->{'_inputFile'};
}

sub percentId{
	my $self=shift;
	$self->{'_percentId'}=shift // return $self->{'_percentId'};
}

sub minimumCharacters{
	my $self=shift;
	$self->{'_minimumDifferentCharacters'}=shift // return $self->{'_minimumDifferentCharacters'};
}

sub maximumMissing{
	my $self=shift;
	$self->{'_maximumMissing'}=shift // return $self->{'_maximumMissing'};
}

sub minimumPresent{
	my $self=shift;
	$self->{'_minimumPresent'}=shift // return $self->{'_minimumPresent'};
}




