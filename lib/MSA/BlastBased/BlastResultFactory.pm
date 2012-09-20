#!/usr/bin/perl
package MSA::BlastBased::BlastResultFactory;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use MSA::BlastBased::BlastResultObject;
use FileInteraction::Fasta::SequenceName;
use feature qw{ switch };

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

sub resultObj{
	my $self=shift;
	$self->{'_resultObj'}=shift // return $self->{'_resultObj'};
}

sub fileHandle{
	my $self=shift;
	$self->{'_fileHandle'}=shift // return $self->{'_fileHandle'};
}

#methods
sub _initialize {
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
	my $first_hsp=0;
	while ( my $line = $self->fileHandle->getline ) {

		if($line =~ m/<\/Hit>/){
			$first_hsp=0;
		}

		if($first_hsp == 1){
			next;
			#we only want the first hsp from every hit
			#skip the rest to avoid doing needless assignment calls
		}

		given($line){

			when(/<\/Iteration>/){
				$first_hsp=0;
				return $self->resultObj();
			}

			when(/<\/Hsp>/){
				$first_hsp=1;
				next;
			}

			when(/<Iteration_iter-num>(\d+)<\/Iteration_iter-num>/){
				$self->resultObj( MSA::BlastBased::BlastResultObject->new() );
				$self->resultObj->iteration_num($1);
			}

			when(/<Iteration_query-def>(.+)<\/Iteration_query-def>/){
				$self->resultObj->query_def($1);
				my $name = FileInteraction::Fasta::SequenceName->new($1);
				$self->resultObj->name( $name->name );
			}

			when(/<Iteration_query-len>(.+)<\/Iteration_query-len>/ ){
				$self->resultObj->query_len($1);
			}

			when(/<Hit_def>(.+)<\/Hit_def>/){
				$currentHit = FileInteraction::Fasta::SequenceName->new($1);
				$self->resultObj->addHit( $currentHit->name );
				$self->resultObj->hitHash->{ $currentHit->name}->addParam( 'hit_def', $1 );
			}

			when(/<Hit_accession>(.+)<\/Hit_accession>/){
				$self->resultObj->hitHash->{$currentHit->name}->addParam( 'hit_accession', $1 );
			}

			when(/<Hsp_query-from>(.+)<\/Hsp_query-from>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_query_from', $1 );
			}

			when(/<Hsp_query-to>(.+)<\/Hsp_query-to>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_query_to', $1 );
			}

			when(/<Hsp_hit-from>(.+)<\/Hsp_hit-from>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hit_from', $1 );
			}

			when(/<Hsp_hit-to>(.+)<\/Hsp_hit-to>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hit_to', $1 );
			}

			when(/<Hsp_identity>(.+)<\/Hsp_identity>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_identity', $1 );
			}

			when(/<Hsp_align-len>(.+)<\/Hsp_align-len>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_align_len', $1 );
			}

			when(/<Hsp_qseq>(.+)<\/Hsp_qseq>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_qseq', $1 );
			}

			when(/<Hsp_hseq>(.+)<\/Hsp_hseq>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hseq', $1 );
			}

			when(/<Hsp_query-frame>(.+)<\/Hsp_query-frame>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_query_frame', $1 );
			}

			when(/<Hsp_hit-frame>(.+)<\/Hsp_hit-frame>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hit_frame', $1 );
			}

			when(/<Hit_len>(.+)<\/Hit_len>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hit_len', $1 );
			}

			when(/<Hsp_score>(.+)<\/Hsp_score>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_score', $1 );
			}

			when(/<Hsp_evalue>(.+)<\/Hsp_evalue>/){
				$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_evalue', $1 );
			}
		}# end of given		
	}    #end of while
}    #end of sub

1;
