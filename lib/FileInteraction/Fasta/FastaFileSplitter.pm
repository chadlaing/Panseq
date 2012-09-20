#!/usr/bin/perl
package FileInteraction::Fasta::FastaFileSplitter;


=pod

=head1 NAME

FileInteraction::Fasta::FastaFileSplitter - Splits a single fasta file into a given number of ~equal sized output fasta files.

=head1 SYNOPSIS



=head1 DESCRIPTION



=head2 Methods


=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing gmail com)

=cut

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Carp;
use IO::File;
use FileInteraction::Fasta::SequenceRetriever;
use Log::Log4perl;
use parent 'FileInteraction::FlexiblePrinter';

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

#class variables
sub _numberOfTempFiles {
	my $self = shift;
	$self->{'_FastaFileSplitter_numberOfTempFiles'} = shift // return $self->{'_FastaFileSplitter_numberOfTempFiles'};
}

sub arrayOfSplitFiles {
	my $self = shift;
	$self->{'_FastaFileSplitter_arrayOfSplitFiles'} = shift // return $self->{'_FastaFileSplitter_arrayOfSplitFiles'};
}

sub logger{
	my $self=shift;
	$self->{'_FastaFileSplitter_logger'} = shift // return $self->{'_FastaFileSplitter_logger'};
}

#methods
sub _initialize {
	my $self = shift;

	#inheritance
	$self->SUPER::_initialize(@_);

	$self->arrayOfSplitFiles( [] );    #init as an anonymous array

	#logging
	$self->logger(Log::Log4perl->get_logger());
}

sub splitFastaFile {
	my $self           = shift;
	my $fileName       = shift;
	my $numberOfSplits = shift;
	
	$self->logger->info("Splitting $fileName into $numberOfSplits files");

	my $fileHandle = IO::File->new( '<' . $fileName ) or die "cannot open $fileName" . $!;
	my %fastaSizes;
	my $currentHeader;
	my $currentSize;

	while ( my $line = $fileHandle->getline ) {
		next if $line eq '';

		if ( $line =~ /^>(.+)/ ) {
			$fastaSizes{$currentHeader} = $currentSize if defined $currentHeader;
			$currentHeader              = $1;
			$currentSize                = 0;
		}
		elsif ( $line =~ /^([\w\-])/ ) {
			my $seqString = $1;
			$currentSize .= length($seqString);
		}
	}
	$fastaSizes{$currentHeader} = $currentSize;

	#avoid creating blank temp files
	my $numberOfSeqs = scalar keys %fastaSizes;
	if ( $numberOfSeqs >= $numberOfSplits ) {
		$self->_numberOfTempFiles($numberOfSplits);
	}
	else {
		$self->_numberOfTempFiles($numberOfSeqs);
	}

	$self->logger->info( "Number of temp query files to create: " . $self->_numberOfTempFiles );

	#package into approximately equal temp files
	#use the Bio::DB in SequenceRetriever
	my $queryDB = FileInteraction::Fasta::SequenceRetriever->new(	
		'inputFile'=>$fileName,
		'outputFile'=>$fileName . '_sequenceSplitterDB'
	);

	my $tempNum    = 0;
	my $count      = 0;
	my @sortedKeys = sort { $fastaSizes{$a} <=> $fastaSizes{$b} } keys %fastaSizes;
	for ( my $start = 0 ; $start < $self->_numberOfTempFiles ; $start++ ) {
		my $tempFileName = $fileName . $start . '.FastaTemp';
		push @{ $self->arrayOfSplitFiles }, $tempFileName;
		$self->outputFH(IO::File->new('>' . $tempFileName or die "Cannot open $tempFileName"));
		
		for(my $seqNum=$start;$seqNum < scalar(@sortedKeys); $seqNum+=$self->_numberOfTempFiles){
			my $seq = $sortedKeys[$seqNum];
			$self->print( '>' . $seq . "\n" . $queryDB->extractRegion($seq) . "\n" );
		}
	}
}

1;
