#!/usr/bin/perl

=pod

=head1 NAME

FileInteraction::SequenceRetriever - Grabs fasta segments from a specified inputFile and returns them or prints them as a multi-fasta file.

=head1 SYNOPSIS

	use FindBin::libs;
	use FileInteraction::SequenceRetriever;
	
	my $retriever = FileInteraction::SequenceRetriever->new(
		'inputFile'=> '/input/file.fasta',
		'outputFile'=> '/output/file.fasta'
	);

=head1 DESCRIPTION

The module allows the easy extraction and printing of a number of sequences and subsequences stored in a hash in the 
format $hashRef->{id}=,3..34,45..199,455..902 etc. It uses Bio::DB::Fasta to create a database of the fasta file specified
as inputFile. Sequences and subsequences can be extracted individually using ->extractRegion or a hash reference can be passed
to extractAndPrintRegionsFromHash to automatically create an output multi-fasta file of the regions of interest.
This module uses Bio::SeqIO and Bio::Seq to manipulate the fasta files and sequences.
This module also inherits from FileInteraction::FlexiblePrinter the methods ->print and ->outputFH. The former prints out to the latter.

=head2 Methods

=head3 inputFile

The absolute location of the fasta file to retrieve sequence from.

=head3 outputFile

The absolute location of the extracted multi-fasta regions file.
If not specified, defaults to STDOUT.

=head3 _setInputFile

Private method to set the read-only variable inputFile upon initialization.

=head3 logger

Stores the logger instance.

=head3 _initialize

Called upon object construction.
Creates logger instance.
Creates a database from the input fasta file if not already created.
Assigns required inputFile and optional outputFile values.

=head3 _createDatabase

Internal method to create a fasta database from inputFile.
Uses Bio::DB::Fasta
Ensures compatability with the Bio::DB::Fasta module.
No lines over 65,536 characters and all lines the same length except the last per sequence.
Uses Bio::SeqIO for this purpose.

=head3 _databaseName

Internal method to store the created databaseName created from the formatted fasta file taken from inputFile.

=head3 extractAndPrintRegionsFromHash

	$obj->extractAndPrintRegionsFromHash(
		'hashRef' 		=> \%hash,
		'cutoffSize'    => 500
	);
	
Takes in a hashRef in the form:
$hashRef->id->,1..5,8..134,678..45999 etc.
The id is the fasta id indexed in the database and the string of locations ,X..Y represent a single subsequence to return.
This method calls the extractRegion method and prints the results greater than or equal to 'cutoffSize' to outputFile.

=head3 extractRegion

	$obj->extractRegion('id'); #form 1
	$obj->extractRegion('id',<start>,<stop>); #form 2

This method returns the entire sequence referenced in the database by 'id' in form 1.
In form 2 it returns the subsequence of 'id' specified by the <start> and <stop> coords.

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing gmail com)

=cut


package FileInteraction::Fasta::SequenceRetriever;

use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use File::Temp;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use Log::Log4perl;

use parent 'FileInteraction::FlexiblePrinter';

sub new{
	my($class)  = shift;
    my $self= $class->SUPER::new(@_); #this calls SequenceRetriever::_initialize, so need to call the SUPER::_initialize
    return $self;
}

sub outputFile{
	my $self=shift;
	$self->{'_outputFile'} = shift // return $self->{'_outputFile'};	
}


sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};	
}

sub _databaseName{
	my $self=shift;
	$self->{'_databaseName'} = shift // return $self->{'_databaseName'};	
}

##methods
sub _initialize{
	my($self)=shift;
    
     #inheritance
    $self->SUPER::_initialize(@_); 

    #logging
    $self->logger(Log::Log4perl->get_logger());      
    
    #outputFile and inputFile should be read only
    my %init=@_; 

    if($init{'outputFile'}){
    	$self->outputFile($init{'outputFile'});
    }
    else{
    	$self->outputFile('');
    	#allows STDOUT to be an option, rather than requiring an outputFile name
    }
    $self->logger->info('outputFile: ' . $self->outputFile);
    $self->_setInputFile($init{'inputFile'}) // confess("Input file required in Sequence Retriever");
    
    #init database
    $self->_createDatabase();

}

sub _setInputFile{
	my $self = shift;
	$self->{'_inputFile'} = shift // confess("_setInputFile requires a file name to be set");
}

sub inputFile{
	my $self=shift;
	return $self->{'_inputFile'};
}

sub _createDatabase{
	my($self)=shift;
	
	if(defined $self->_databaseName){
		#database already created, no need to duplicte
		$self->logger->info("Database already created" );
	}
	else{
		$self->logger->info("Creating database");
		
		my $originalFH = Bio::SeqIO->new(-file=>'<'.$self->inputFile, -format=>'fasta') or die "$!";
		
		my $dbFileName = $self->outputFile . '_dbtemp';
		
		$self->logger->debug($dbFileName);
		my $outputFH = Bio::SeqIO->new(-file=>'>'. $dbFileName, -format=>'fasta') or die "$!";
		$outputFH->width(80);
		
		while(my $seq = $originalFH->next_seq()){
			$outputFH->write_seq($seq);
		}
		$originalFH->close();
		$outputFH->close();
		
		$self->_databaseName(Bio::DB::Fasta->new($dbFileName));
	}	
}

sub extractAndPrintRegionsFromHash{
	my($self)=shift;
	
	my %params=@_;
	
	my $hashRef = $params{'hashRef'} // confess ('Hash reference required in extractAndPrintRegionsFromHash');
	my $cutoffSize = $params{'cutoffSize'} // confess ('cutoffSize required in extractAndPrintRegionsFromHash');
	
	#expects hashvalues of $hash->id->,1..5,8..134,678..45999 etc. 
	#same format as generated from NovelRegionFinder
	
	#open fileHandle for printing if specified; otherwise FlexiblePrinter defaults to STDOUT
	if($self->outputFile ne ''){
		$self->outputFH(
    		(IO::File->new('>'.$self->outputFile) or die "Could not open ". $self->outputFile . " for printing. $!\n")
    	);
	}
		
	foreach my $id(keys %{$hashRef}){
		my $coordString = $hashRef->{$id};
			
		my $novelCounter=0;
		while($coordString =~ /\,(\d+)\.\.(\d+)/gc){				
			my $start =$1;
			my $end =$2;
			my $length=$end-$start+1;
				
			if(defined $cutoffSize){
				next unless $length >= $cutoffSize;
			}
				
			$novelCounter++;
				
			$self->print('>' . $id . '|NovelRegion=' . $novelCounter . '|Start=' . $start . '|End=' . $end . '|Length=' . $length . "\n" . $self->extractRegion($id,$start,$end) . "\n");	
		}# end while
	}# end foreach
	$self->outputFH->close();	
}

sub extractRegion{
	my($self)=shift;
	
	if(@_){
		my $seqID=shift;
		my $startBp;
		my $endBp;
		
		if(scalar(@_)==2){
			$startBp=shift;
			$endBp=shift;
		}

		my $seqObj = $self->_databaseName->get_Seq_by_id($seqID) // confess ("Cannot locate $seqID");
		$self->logger->debug("Found $seqID\n");
		if(defined $startBp && defined $endBp){
			return $seqObj->subseq($startBp=>$endBp);
		}
		else{
			return $seqObj->seq;
		}					
	}
	else{
		print STDERR "extract region is missing parameters!\n";
		exit(1);
	}
}

1;
