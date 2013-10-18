#!/usr/bin/env perl
package Modules::Alignment::BlastResults;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Log::Log4perl;
use Carp;
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


sub _outFH{
	my $self=shift;
	$self->{'__outFH'} = shift // return $self->{'__outFH'};
}

sub _storedQid{
	my $self=shift;
	$self->{'__storedQid'} = shift // return $self->{'__storedQid'};
}

sub _storedResults{
	my $self=shift;
	$self->{'__storedResults'} = shift // return $self->{'__storedResults'};
}



#methods
sub _initialize {
	my ($self) = shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::BlastResults");
	my $outFile = shift // $self->logger->logdie("No file sent to BlastResults");
	my $cutoff = shift // $self->logger->logdie("percentIdentityCutoff required");

	$self->percentIdentityCutoff($cutoff);
	$self->_outFH(IO::File->new('<' . $outFile)) // $self->logger->logdie("$!");
	$self->_storedResults({});
	$self->_storedQid(undef);
}

sub percentIdentityCutoff{
	my $self=shift;
	$self->{'_percentIdentityCutoff'} = shift // return $self->{'_percentIdentityCutoff'};
}

sub getNextResult{
	my $self=shift;
		
	my $results;
	my $line = $self->_outFH->getline();
	while($line){
		my $nextLine = $self->_outFH->getline();	
		
		if(!defined $nextLine){
			return $results;
		}
		
		$line =~ s/\R//g;
		$nextLine =~ s/\R//g;
		my @la = split("\t",$line);
		my @nextLa = split("\t",$nextLine);
		
		if($la[1] ne $nextLa[1]){
			return $results;
		}

		#'outfmt'=>'"6 
		# [0]sseqid 
		# [1]qseqid 
		# [2]sstart 
		# [3]send 
		# [4]qstart 
		# [5]qend 
		# [6]slen 
		# [7]qlen 
		# [8]pident 
		# [9]length"',
		# [10]sseq,
		# [11]qseq

		my $sName = Modules::Fasta::SequenceName->new($la[0]);
		my $sNameName= $sName->name;			

		$self->_getPercentId($la[7],$la[8],$la[4],$la[5]);			
			
		
	}	
}

sub _getPercentId{
	my $self=shift;
	my $qlen =shift // $self->logger->logdie("Missing qlen");
	my $pident =shift // $self->logger->logdie("Missing pident");
	my $qstart = shift // $self->logger->logdie("Missing qstart");
	my $qend =shift // $self->logger->logdie("Missing qend");

	my $queryHitLength = abs($qend - $qstart) +1;
	
	my $percentId = $queryHitLength / $qlen * $pident;
	 #$self->logger->info("queryHitLength: $queryHitLength");
	 # $self->logger->info("qlen: $qlen");
	 #$self->logger->info("pident: $pident");
	 #$self->logger->info("percentId: $percentId");
	return $percentId;
}

1;
