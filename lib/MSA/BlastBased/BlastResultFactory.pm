#!/usr/bin/perl
package BlastResultFactory;

use FindBin::libs;
use IO::File;
use MSA::BlastBased::BlastResultObject;
use FileInteraction::Fasta::SequenceName;

use Object::Tiny::RW qw{
  resultObj
  fileHandle
};

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->initialize(@_);
	return $self;
}

#methods
sub initialize {
	my ($self) = shift;

	if (@_) {
		my $fileHandle = shift;
		$self->fileHandle($fileHandle);
	}
	else {
		print STDERR "no filename specified!\n";
		exit(1);
	}
}

sub nextResult {
	my ($self) = shift;

	my $currentHit;
	while ( my $line = $self->fileHandle->getline ) {
		if ( $line =~ m/<Iteration_iter-num>(\d+)<\/Iteration_iter-num>/ ) {
			$self->resultObj( BlastResultObject->new() );
			$self->resultObj->iteration_num($1);
		}

		if ( $line =~ m/<Iteration_query-def>(.+)<\/Iteration_query-def>/ ) {
			$self->resultObj->query_def($1);
			my $name = SequenceName->new($1);
			$self->resultObj->name( $name->name );
		}

		if ( $line =~ m/<Iteration_query-len>(.+)<\/Iteration_query-len>/ ) {
			$self->resultObj->query_len($1);
		}

		if ( $line =~ m/<Hit_def>(.+)<\/Hit_def>/ ) {
			$currentHit = SequenceName->new($1);
			$self->resultObj->addHit( $currentHit->name );
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hit_def', $1 );
		}

		if ( $line =~ m/<Hsp_query-from>(.+)<\/Hsp_query-from>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_query_from', $1 );
		}

		if ( $line =~ m/<Hsp_query-to>(.+)<\/Hsp_query-to>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_query_to', $1 );
		}

		if ( $line =~ m/<Hsp_hit-from>(.+)<\/Hsp_hit-from>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_hit_from', $1 );
		}

		if ( $line =~ m/<Hsp_hit-to>(.+)<\/Hsp_hit-to>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_hit_to', $1 );
		}

		if ( $line =~ m/<Hsp_identity>(.+)<\/Hsp_identity>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_identity', $1 );
		}

		if ( $line =~ m/<Hsp_align-len>(.+)<\/Hsp_align-len>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_align_len', $1 );
		}

		if ( $line =~ m/<Hsp_qseq>(.+)<\/Hsp_qseq>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_qseq', $1 );
		}

		if ( $line =~ m/<Hsp_hseq>(.+)<\/Hsp_hseq>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_hseq', $1 );
		}

		if ( $line =~ m/<Hsp_query-frame>(.+)<\/Hsp_query-frame>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_query_frame', $1 );
		}

		if ( $line =~ m/<Hsp_hit-frame>(.+)<\/Hsp_hit-frame>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hsp_hit_frame', $1 );
		}

		if ( $line =~ m/<Hit_len>(.+)<\/Hit_len>/ ) {
			$self->resultObj->hitHash->{ $currentHit->name }
			  ->addParam( 'hit_len', $1 );
		}

		if ( $line =~ m/<\/Iteration>/ ) {
			return $self->resultObj();
		}
	}    #end of while
}    #end of sub

1;
