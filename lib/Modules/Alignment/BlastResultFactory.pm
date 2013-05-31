#!/usr/bin/env perl
package Modules::Alignment::BlastResultFactory;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Log::Log4perl;
use Carp;
use Modules::Alignment::BlastResult;
use Modules::Fasta::SequenceName;

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
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

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::BlastResultFactory\n");

	my $fileHandle = shift // $self->logger->logdie("Modules::Alignment::BlastResultFactory requires a file handle");
	$self->fileHandle($fileHandle);
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

		if($line =~ m/<\/Iteration>/){
			$first_hsp=0;
			return $self->resultObj();
		}
		elsif($line =~ m/<\/Hsp>/){
			$first_hsp=1;
			next;
		}
		elsif($line =~ m/<Iteration_iter-num>(\d+)<\/Iteration_iter-num>/){
			$self->resultObj( Modules::Alignment::BlastResult->new() );
			$self->resultObj->iteration_num($1);
		}
		elsif($line =~ m/<Iteration_query-def>(.+)<\/Iteration_query-def>/){
			$self->resultObj->query_def($1);
			my $name = Modules::Fasta::SequenceName->new($1);
			$self->resultObj->name( $name->name );
		}
		elsif($line =~ m/<Iteration_query-len>(.+)<\/Iteration_query-len>/ ){
			$self->resultObj->query_len($1);
		}
		elsif($line =~ m/<Hit_def>(.+)<\/Hit_def>/){
			$currentHit = Modules::Fasta::SequenceName->new($1);
			$self->resultObj->addHit( $currentHit->name );
			$self->resultObj->hitHash->{ $currentHit->name}->addParam( 'hit_def', $1 );
		}
		elsif($line =~ m/<Hit_accession>(.+)<\/Hit_accession>/){
			$self->resultObj->hitHash->{$currentHit->name}->addParam( 'hit_accession', $1 );
		}
		elsif($line =~ m/<Hsp_query-from>(.+)<\/Hsp_query-from>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_query_from', $1 );
		}
		elsif($line =~ m/<Hsp_query-to>(.+)<\/Hsp_query-to>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_query_to', $1 );
		}
		elsif($line =~ m/<Hsp_hit-from>(.+)<\/Hsp_hit-from>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hit_from', $1 );
		}
		elsif($line =~ m/<Hsp_hit-to>(.+)<\/Hsp_hit-to>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hit_to', $1 );
		}
		elsif($line =~ m/<Hsp_identity>(.+)<\/Hsp_identity>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_identity', $1 );
		}
		elsif($line =~ m/<Hsp_align-len>(.+)<\/Hsp_align-len>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_align_len', $1 );
		}
		elsif($line =~ m/<Hsp_qseq>(.+)<\/Hsp_qseq>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_qseq', $1 );
		}
		elsif($line =~ m/<Hsp_hseq>(.+)<\/Hsp_hseq>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hseq', $1 );
		}
		elsif($line =~ m/<Hsp_query-frame>(.+)<\/Hsp_query-frame>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_query_frame', $1 );
		}
		elsif($line =~ m/<Hsp_hit-frame>(.+)<\/Hsp_hit-frame>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_hit_frame', $1 );
		}
		elsif($line =~ m/<Hit_len>(.+)<\/Hit_len>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hit_len', $1 );
		}
		elsif($line =~ m/<Hsp_score>(.+)<\/Hsp_score>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_score', $1 );
		}
		elsif($line =~ m/<Hsp_evalue>(.+)<\/Hsp_evalue>/){
			$self->resultObj->hitHash->{ $currentHit->name }->addParam( 'hsp_evalue', $1 );
		}
	}    #end of while
}    #end of sub

1;
