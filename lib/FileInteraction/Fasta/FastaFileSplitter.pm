#!/usr/bin/perl
package FastaFileSplitter;

use strict;
use warnings;
use FindBin::libs;
use Carp;
use IO::File;
use FileInteraction::Fasta::SequenceRetriever;
use FileInteraction::FileManipulation;
use Logging::Logger;
our @ISA = qw/Logger/;

#object creation
sub new{
	my($class)=shift;
    my $self = {};
    bless ($self, $class);
    $self->_fastaFileSplitterInitialize(@_);
    return $self;
}

#class variables
sub _numberOfTempFiles{
	my $self=shift;
	$self->{'_FastaFileSplitter_numberOfTempFiles'}=shift // return $self->{'_FastaFileSplitter_numberOfTempFiles'};
}

sub arrayOfSplitFiles{
	my $self=shift;
	$self->{'_FastaFileSplitter_arrayOfSplitFiles'}=shift // return $self->{'_FastaFileSplitter_arrayOfSplitFiles'};
}

#methods
sub _fastaFileSplitterInitialize{
	my $self=shift;
	
	$self->arrayOfSplitFiles([]); #init as an anonymous array
	
	#inheritance
	$self->_loggerInitialize();	
}

sub splitFastaFile{
	my $self=shift;
	my $fileName=shift;
	my $numberOfSplits=shift;
	
	$self->logger->info("INFO:\tSplitting $fileName into $numberOfSplits files");	
	
	my $fileHandle = IO::File->new( '<' . $fileName ) or die "cannot open $fileName" . $!;
	my %fastaSizes;
	my $currentHeader;
	my $currentSize;
	
	while(my $line = $fileHandle->getline){
		next if $line eq '';
		
		if($line =~ /^>(.+)/){
			$fastaSizes{$currentHeader}=$currentSize if defined $currentHeader;
			$currentHeader=$1;
			$currentSize=0;
		}
		elsif($line =~ /^([\w\-])/){
			my $seqString = $1;
			$currentSize .= length($seqString);
		}	
	}
	$fastaSizes{$currentHeader}=$currentSize;
	
	#avoid creating blank temp files
	my $numberOfSeqs = scalar keys %fastaSizes;
	if($numberOfSeqs >= $numberOfSplits){
		$self->_numberOfTempFiles($numberOfSplits);
	}
	else{
		$self->_numberOfTempFiles($numberOfSeqs);
	}
	
	$self->logger->info("INFO:\tNumber of temp query files to create: " . $self->_numberOfTempFiles);
	
	#package into approximately equal temp files
	#use the Bio::DB in SequenceRetriever
	my $queryDB = SequenceRetriever->new($fileName);
	
	my $tempNum=0;
	my $count=0;
	my @sortedKeys = sort {$fastaSizes{$a} <=> $fastaSizes{$b}} keys %fastaSizes;
	
	for(my $start=0;$start < $self->_numberOfTempFiles;$start++){
		my $tempFileName = $fileName . $start. '.FastaTemp';
		push @{$self->arrayOfSplitFiles}, $tempFileName;
		my $tempFH = IO::File->new('>' . $tempFileName) or die "Cannot open $tempFileName\n $!";
		
		for(my $seqNum=$start;($seqNum + $self->_numberOfTempFiles) < scalar(@sortedKeys); $seqNum+=$self->_numberOfTempFiles){
			my $seq=$sortedKeys[$seqNum];
			print $tempFH ('>' . $seq . "\n" . $queryDB->extractRegion($seq) . "\n");
		
			$self->logger->debug("DEBUG: Adding $seq to $tempFileName");			
		}
		$tempFH->close();		
	}
}

1;