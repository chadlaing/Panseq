#!/usr/bin/perl

package SegmentMaker;

use FindBin::libs;
use IO::File;
use Logging::Logger;
use FileInteraction::FlexiblePrinter;

our @ISA = qw{FlexiblePrinter Logger}; #includes outputFilehandle and printOut

use Object::Tiny::RW;

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_segmentMakerInitialize(@_);
    return $self;
}

##methods
sub _segmentMakerInitialize{
	my($self)=shift;
	
	if(@_){
		my $fh=shift;
		$self->outputFilehandle($fh);
	}
	
	#inheritance
	$self->_loggerInitialize();	
}

sub segmentTheSequence{
	my($self)=shift;
	
	if(scalar(@_)==2){
		my $inputFile=shift;
		my $fragmentSize = shift; 
		
		$self->logger->info("INFO\tSegmenting the file: $inputFile into $fragmentSize bp fragments");
		
		open(INPUTSEG, "<", $inputFile) or die "no input file $!";		

		local $/=">";
		my $fastaHeader;
		my $fastaSequence;
		
		while(my $fastaSegment=<INPUTSEG>){
			
			next if $fastaSegment eq '>';
			
			#get header and sequence
			if($fastaSegment =~ m/^(.*)\n/){
				$fastaHeader=$1;
				#$fastaHeader =~ s/\s/_/g;
				$fastaSegment =~ s/^.+\n//;
				$fastaSegment =~ s/[\W>]//g;
				$fastaSequence = $fastaSegment;
				$fastaSegment='';
			}
			else{
				print STDERR "not a fasta segment in segmentSequence!\n";
				exit(1);
			}
			
			#segment the sequence
			my $endPosition=0;
			my $startPosition=0;
			my $segmentCounter=0;
			my $currentStart=0;
			my $sequenceLength = length($fastaSequence);
			while($sequenceLength > $currentStart){
				$segmentCounter++;
				my $sequenceToPrint;
				#if twice fragmentation size, print fragmentation size, else print whole thing
				if($sequenceLength >= (2*$fragmentSize) + $currentStart){				
					$sequenceToPrint = substr($fastaSequence, $currentStart, $fragmentSize);	
					$currentStart = $currentStart + $fragmentSize;			
				}
				else{
					$sequenceToPrint = substr($fastaSequence, $currentStart);
					$currentStart = $sequenceLength +1;
				}
				
				#update start/end values
				$startPosition=$endPosition+1;
				$endPosition= $endPosition + length($sequenceToPrint);
				
				#print em out
				$self->printOut('>' . $fastaHeader . '|Segment=' . $segmentCounter . '|SegmentLength=' . length($sequenceToPrint) . '|Start=' . $startPosition . '|End=' . $endPosition . "\n");
				$self->printOut($sequenceToPrint . "\n");
				my $sequenceLength = length($fastaSequence);
		}
	}
	
	close INPUTSEG;
	}
	else{
		print STDERR "segmentation parameters not defined!";
		exit(1);
	}
}

1;
