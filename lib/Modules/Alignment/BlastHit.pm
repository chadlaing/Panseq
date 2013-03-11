#!/usr/bin/env perl
package Modules::Alignment::BlastHit;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use Log::Log4perl;
use Carp;

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
	$self->logger->debug("Logger initialized in Modules::Alignment::BlastHit\n");
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

sub hit_def{
	my $self=shift;
	$self->{'_hit_def'}=shift // return $self->{'_hit_def'};
}

sub hit_accession{
	my $self=shift;
	$self->{'_hit_accession'}=shift // return $self->{'_hit_accession'};
}

sub hit_len{
	my $self=shift;
	$self->{'_hit_len'}=shift // return $self->{'_hit_len'};
}

sub hsp_query_from{
	my $self=shift;
	$self->{'_hsp_query_from'}=shift // return $self->{'_hsp_query_from'};
}

sub hsp_query_to{
	my $self=shift;
	$self->{'_hsp_query_to'}=shift // return $self->{'_hsp_query_to'};
}

sub hsp_hit_from{
	my $self=shift;
	$self->{'_hsp_hit_from'}=shift // return $self->{'_hsp_hit_from'};
}

sub hsp_hit_to{
	my $self=shift;
	$self->{'_hsp_hit_to'}=shift // return $self->{'_hsp_hit_to'};
}

sub hsp_identity{
	my $self=shift;
	$self->{'_hsp_identity'}=shift // return $self->{'_hsp_identity'};
}

sub hsp_align_len{
	my $self=shift;
	$self->{'_hsp_align_len'}=shift // return $self->{'_hsp_align_len'};
}

sub hsp_qseq{
	my $self=shift;
	$self->{'_hsp_qseq'}=shift // return $self->{'_hsp_qseq'};
}
sub hsp_hseq{
	my $self=shift;
	$self->{'_hsp_hseq'}=shift // return $self->{'_hsp_hseq'};
}

sub hsp_query_frame{
	my $self=shift;
	$self->{'_hsp_query_frame'}=shift // return $self->{'_hsp_query_frame'};
}

sub hsp_hit_frame{
	my $self=shift;
	$self->{'_hsp_hit_frame'}=shift // return $self->{'_hsp_hit_frame'};
}

sub hsp_score{
	my $self=shift;
	$self->{'_hsp_score'}=shift // return $self->{'_hsp_score'};
}

sub hsp_evalue{
	my $self=shift;
	$self->{'_hsp_evalue'}=shift // return $self->{'_hsp_evalue'};
}

#methods
sub addParam {
	my ($self) = shift;

	if ( scalar(@_) == 2 ) {
		my $param = shift;
		my $value = shift;

		unless ( defined $self->$param ) {
			$self->$param($value);
		}
	}
	else {
		$self->logger->logconfess("missing values to add param!");
	}
}

sub setSequence {
	my ($self) = shift;

	if ( defined $self->hsp_query_frame && defined $self->hsp_hit_frame ) {
		if ( $self->hsp_query_frame ne $self->hsp_hit_frame ) {
			$self->switchHitToFrom();
		}
	}
	else {
		$self->logger->logconfess("query frames not defined");
	}
}

sub switchHitToFrom {
	my ($self) = shift;

	my $sQ = $self->hsp_hit_from;
	my $eQ = $self->hsp_hit_to;
	$self->hsp_hit_from($eQ);
	$self->hsp_hit_to($sQ);
}

1;
