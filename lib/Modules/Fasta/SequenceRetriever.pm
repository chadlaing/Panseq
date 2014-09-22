#!/usr/bin/perl

=pod

=head1 NAME

Modules::Fasta::SequenceRetriever - Grabs fasta segments from a specified inputFile and returns them or prints them as a multi-fasta file.

=head1 SYNOPSIS

	use FindBin::libs;
	use Modules::Fasta::SequenceRetriever;
	
	my $retriever = Modules::Fasta::SequenceRetriever->new(
		'inputFile'=> '/input/file.fasta',
		'databaseFile'=> 'queryfile_dbtemp'
	);

=head1 DESCRIPTION

The module allows the easy extraction and printing of a number of sequences and subsequences stored in a hash in the 
format $hashRef->{id}=,3..34,45..199,455..902 etc. It uses Bio::DB::Fasta to create a database of the fasta file specified
as inputFile. Sequences and subsequences can be extracted individually using ->extractRegion or a hash reference can be passed
to extractAndPrintRegionsFromHash to automatically create an output multi-fasta file of the regions of interest.
This module uses Bio::SeqIO and Bio::Seq to manipulate the fasta files and sequences.
This module also uses the role Roles::FlexiblePrinter the methods ->printOut and ->outputFH. The former prints out to the latter.

=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing gmail com)

=cut

package Modules::Fasta::SequenceRetriever;

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

sub new {
	my $class = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


=head3 _initialize

Called upon object construction.
Creates logger instance.
Creates a database from the input fasta file if not already created.
Assigns required inputFile and databaseFile value.

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 
    $self->logger->debug("Logger initialized in Modules::Fasta::SequenceRetriever");     
    
    #inputFile should be read only
    my %params=@_; 

    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::Fasta::SequenceRetriever");
		}
	}	

    #init database
    $self->_createDatabase();
}

=head2 databaseFile

Specify the file location for creation of the database.

=cut

sub databaseFile{
	my $self=shift;
	$self->{'_databaseFile'}=shift // return $self->{'_databaseFile'};
}

=head3 logger

Stores the logger instance.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};	
}


=head3 _databaseName

Internal method to store the created databaseName created from the formatted fasta file taken from inputFile.

=cut

sub _databaseName{
	my $self=shift;
	$self->{'_databaseName'} = shift // return $self->{'_databaseName'};	
}



=head3 inputFile

The absolute location of the fasta file to retrieve sequence from.

=cut

sub inputFile{
	my $self=shift;
	$self->{'_inputFile'}=shift // return $self->{'_inputFile'};
}

=head3 cutoffSize

The minimum size a retrieved sequence must be to avoid being discarded.
If left undefined, there is no limit.

=cut

sub cutoffSize{
	my $self=shift;
	$self->{'_cutoffSize'} = shift // return $self->{'_cutoffSize'};
}



=head3 _createDatabase

Internal method to create a fasta database from inputFile.
Uses Bio::DB::Fasta
Ensures compatability with the Bio::DB::Fasta module.
No lines over 65,536 characters and all lines the same length except the last per sequence.
Uses Bio::SeqIO for this purpose.

=cut

sub _createDatabase{

	my($self)=shift;
	
	if(defined $self->_databaseName){
		#database already created, no need to duplicte
		$self->logger->info("Database already created" );
	}
	else{
		$self->logger->info("Creating database");
		
		my $originalFH = Bio::SeqIO->new(-file=>'<'.$self->inputFile, -format=>'fasta') or die "$!";
		my $dbFileName = $self->databaseFile;

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



=head3 extractRegion

	$obj->extractRegion('id'); #form 1
	$obj->extractRegion('id',<start>,<stop>); #form 2

This method returns the entire sequence referenced in the database by 'id' in form 1.
In form 2 it returns the subsequence of 'id' specified by the <start> and <stop> coords.

=cut

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

		my $seqObj = $self->_databaseName->get_Seq_by_id($seqID) // $self->logger->logconfess("Cannot locate $seqID in database based on " . $self->databaseFile);
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
