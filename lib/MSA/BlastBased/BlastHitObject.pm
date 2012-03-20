#!/usr/bin/perl
package BlastHitObject;

use FindBin::libs;
use FileInteraction::Fasta::SequenceName;

use Object::Tiny::RW qw{
  hit_def
  hit_len
  hsp_query_from
  hsp_query_to
  hsp_hit_from
  hsp_hit_to
  hsp_identity
  hsp_align_len
  hsp_qseq
  hsp_hseq
  hsp_query_frame
  hsp_hit_frame
};

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
		print STDERR "missing values to add param!\n";
		exit(1);
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
		print STDERR "query frames not defined\n";
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
