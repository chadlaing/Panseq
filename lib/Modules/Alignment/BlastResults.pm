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


#we want to store UNDEF in this variable, which makes the standard way always return the stored value, rather than setting it to UNDEF.
sub _setStoredLine{
	my $self=shift;
	$self->{'__storedLine'} = shift;	
}

sub _getStoredLine{
	my $self = shift;
	return $self->{'__storedLine'}
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
	
	#set the first line in _storedLine to prevent infinite looping
	$self->_setStoredLine($self->_outFH->getline());
}

sub percentIdentityCutoff{
	my $self=shift;
	$self->{'_percentIdentityCutoff'} = shift // return $self->{'_percentIdentityCutoff'};
}

sub getNextResult{
	my $self=shift;
		
	my $results;
	my $line = $self->_getStoredLine(); #on object creation, this is initialized with the first line
	$self->logger->info("Calling getNextResult with line: $line");
	my $counter=0;
	while($line){		
		$counter++;
		my $nextLine = $self->_outFH->getline();
		$self->_setStoredLine($nextLine);	
		
		my @nextLa;
		if(defined $nextLine){
			$nextLine =~ s/\R//g;
			@nextLa = split("\t",$nextLine);
		}
		else{
			$self->logger->info("nextLine UNDEF");
		}
		
		$line =~ s/\R//g;		
		my @la = split("\t",$line);		
		
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
		
		if($self->_getPercentId($la[7],$la[8],$la[4],$la[5]) > $self->percentIdentityCutoff){
			if(defined $results && defined $results->{$sNameName}){
				#nothing
			}
			else{
				$self->logger->debug("Adding $sNameName to results");
				$results->{$sNameName}=\@la;
			}				
		}		
		
		if(!defined $nextLa[0] || $la[1] ne $nextLa[1]){
			$self->logger->info("nextLa not defined or la1 $la[1] and nextLa1 $nextLa[1] not equal");
			return $results;
		}
		$line = $nextLine;
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
