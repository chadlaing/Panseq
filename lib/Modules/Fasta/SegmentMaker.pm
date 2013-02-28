#!/usr/bin/env perl

=pod

=head1 NAME

FileInteraction::Fasta::SegmentMaker

=head1 SYNOPSIS

	use FindBin::libs;
	use FileInteraction::Fasta::SegmentMaker;
	
	my $segmenter = SegmentMaker->new(
		'outputFile' => '/my/outputfile.txt',
		'inputFile'  => '/my/inputfile.txt',
		'segmentSize' => 500
	);
	
	$segmenter->segmentTheSequence();

=head1 DESCRIPTION

This module takes a fasta file, fragments it into user-defined segments and outputs these as a single 
multi-fasta file. In cases where the size of a fasta sequence is not evenly divisible by the segment
size, then the remainder is added to the last fragment, such that all segments will be of N size, with
the last segment being between N and 2N-1.

This module uses Bio::SeqIO and Bio::Seq to manipulate the fasta files and sequences.

=head2 Methods

=head3 inputFile

The absolute location of the fasta file being segmented.

=head3 outputFile

The absolute location of the segmented fasta file.

=head3 logger

Stores the logger instance.

=head3 segmentSize

Sets the size that inputFile will be segmented into.

=head3 segmentTheSequence

The function that actually does the work.

=head3 _initialize

Called upon object construction.
Initializes inputFile, outputFile, segmentSize and creates a logger instance.

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing gmail com)

=cut


package Modules::Fasta::SegmentMaker;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Bio::SeqIO;
use Bio::Seq;
use Log::Log4perl;
use Carp;

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_initialize(@_);
    return $self;
}

#variables
sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub inputFile{
	my $self=shift;
	$self->{'__inputFH'} = shift // return $self->{'__inputFH'};	
}

sub outputFile{
	my $self=shift;
	$self->{'__outputFile'} = shift // return $self->{'__outputFile'};	
}

sub segmentSize{
	my $self=shift;
	$self->{'__segmentSize'} = shift // return $self->{'__segmentSize'};	
}


##methods
sub _initialize{
	my($self)=shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());

	my %params=@_;
	
	#on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::Fasta::SegmentMaker");
		}
	}	

	unless(defined $self->segmentSize && defined $self->inputFile && defined $self->outputFile){
		$self->logger->logconfess("Modules::Fasta::SegmentMaker requires:\n 
			\t'segmentSize'\n
			\t'inputFile'\n
			\t'outputFile'\n
			to be defined. Missing one or more parameters"
		);
	}
	
}

sub segmentTheSequence{
	my $self =shift;
		
	$self->logger->info('Segmenting ' . $self->inputFile . ' into ' . $self->segmentSize . 'bp segments');
	my $inFH =Bio::SeqIO->new(-file=>'<'. $self->inputFile, -format=>'fasta');
	my $outputFH = Bio::SeqIO->new(-file=>'>'. $self->outputFile, -format=>'fasta');
	$outputFH->width(80);
		
	while(my $seq = $inFH->next_seq()){
		my $segmentCounter=1;
		#if twice fragmentation size, print fragmentation size, else print whole thing
		for(my $i=1; $i<=$seq->length; $i+=$self->segmentSize){
			my $endCalc = ($i + (2 * $self->segmentSize) -1);
			my $end =  ($endCalc > $seq->length) ? $seq->length : ($endCalc - $self->segmentSize);
			
			my $newId = $seq->id() . $seq->desc();
			$newId =~ s/\s//g;
			$newId .= ('|Segment=' . $segmentCounter . '|SegmentLength=' . ($end -$i +1) . '|Start=' . $i . '|End=' . $end);

			my $tempSeq = Bio::Seq->new(
				-seq => $seq->subseq($i,$end),
				-id => $newId,
				-accession_number => $seq->accession_number,
				#-desc => $newDesc
			);
			$outputFH->write_seq($tempSeq);
			$segmentCounter++;
			last if $end == $seq->length; #this avoids duplicating the last segment of sequence
		}			
	}
	$inFH->close();
	$outputFH->close();
}

1;
