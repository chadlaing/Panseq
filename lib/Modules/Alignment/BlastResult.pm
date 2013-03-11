#!/usr/bin/env perl
package Modules::Alignment::BlastResult;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Log::Log4perl;
use Carp;
use Modules::Alignment::BlastHit;
use Modules::Fasta::SequenceName;

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

sub _initialize{
	my $self=shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::BlastResult\n");


	#init hash
	$self->hitHash({});
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

sub name{
	my $self=shift;
	$self->{'_name'}=shift // return $self->{'_name'};
}

sub firstHitName{
	my $self=shift;
	$self->{'_firstHitName'}=shift // return $self->{'_firstHitName'};
}

sub iteration_num{
	my $self=shift;
	$self->{'_iteration_num'}=shift // return $self->{'_iteration_num'};
}

sub query_def{
	my $self=shift;
	$self->{'_query_def'}=shift // return $self->{'_query_def'};
}

sub query_len{
	my $self=shift;
	$self->{'_query_len'}=shift // return $self->{'_query_len'};
}

sub hitHash{
	my $self=shift;
	$self->{'_hitHash'}=shift // return $self->{'_hitHash'};
}

#methods

sub addHit {
	my ($self) = shift;

	if (@_) {
		my $hitName = shift;

		unless ( defined $self->hitHash && defined $self->hitHash->{$hitName} )
		{
			my $hit  = Modules::Alignment::BlastHit->new();
			my $name = Modules::Fasta::SequenceName->new($hitName);
			$self->firstHitName( $name->name )
			  unless defined $self->firstHitName
			;    #used when any hit object is needed, not a specific one
			#$self->addToHash( 'hitHash', $name->name, $hit );
			$self->hitHash->{$name->name}=$hit;
		}
	}
	else {
		print STDERR "no hit sent!\n";
		exit(1);
	}
}

1;
