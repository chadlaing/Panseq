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
use Role::Tiny::With;

with 'Roles::GetNewIdStartEnd';

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

sub _fragmentHash{
	my $self=shift;
	$self->{'__fragmentHash'} = shift // return $self->{'__fragmentHash'};	
}




##methods
sub _initialize{
	my($self)=shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Fasta::SegmentMaker");

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

	$self->_fragmentHash({});
}

sub segmentTheSequence{
	my $self =shift;
		
	$self->logger->info('Segmenting ' . $self->inputFile . ' into ' . $self->segmentSize . 'bp segments');
	my $inFH =Bio::SeqIO->new(-file=>'<'. $self->inputFile, -format=>'fasta');
	my $outputFH = Bio::SeqIO->new(-file=>'>'. $self->outputFile, -format=>'fasta');
	$outputFH->width(80);
		
	while(my $seq = $inFH->next_seq()){
		#if twice fragmentation size, print fragmentation size, else print whole thing
		#store seq details in hashRef for use in the next_fragment iterator
		$self->_fragmentHash->{sequence}=$seq;
		$self->_fragmentHash->{position}=1;
		my $id = $seq->id() . $seq->desc();
		$id =~ s/\s/_/g;
		$self->_fragmentHash->{newId} = $id;

		while(my $frag = $self->next_fragment()){
			$outputFH->write_seq($frag);
		}		
	}
	$inFH->close();
	$outputFH->close();
}
			

sub next_fragment{
	my $self = shift;

	my $bpRemaining = $self->_fragmentHash->{sequence}->length() - $self->_fragmentHash->{position} + 1;
	
	my $bpToTake;
	if($bpRemaining >= (2 * $self->segmentSize)){
		$bpToTake = $self->segmentSize;
	}
	else{
		$bpToTake = $bpRemaining;
	}

	my $start = $self->_fragmentHash->{position};
	my $end = $start + $bpToTake -1;

	if($start > $end){
		return undef;
	}

	#uses Roles::GetNewIdStartEnd to implement _getNewIdStartEnd
	my($relId, $relStart, $relEnd) = $self->_getNewIdStartEnd($self->_fragmentHash->{newId},$start,$end);
	$relId .= '_(' . $relStart . '..' . $relEnd . ')';

	my $tempSeq = Bio::Seq->new(
		-seq => $self->_fragmentHash->{sequence}->subseq($start,$end),
		-id => $relId,
		-accession_number => $self->_fragmentHash->{sequence}->accession_number,
	);

	$self->_fragmentHash->{position}=$end+1;
	return $tempSeq;
}

1;
