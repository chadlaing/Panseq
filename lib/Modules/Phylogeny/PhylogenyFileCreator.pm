#!/usr/bin/env perl

package Modules::Phylogeny::PhylogenyFileCreator;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use IO::File;
use Log::Log4perl;

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


=head3 _initialize

Initializes the logger.
Assigns all values to class variables.

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 

    $self->logger->info("Logger initialized in Modules::PanGenome::PanGenome");  

    my %params = @_;

    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::PanGenome::PanGenome");
		}
	}	


	#deafults
	$self->_setDefaults();
}


=head2 _setDefaults

Sets the defaults for a run if not specified.

=cut


sub _setDefaults{
	my $self=shift;

	unless(defined $self->outputFormat){
		$self->outputFormat('phylip');
	}


}


=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

sub nameOrderArray{
	my $self=shift;
	$self->{'_nameOrderArray'}=shift // return $self->{'_nameOrderArray'};
}

sub outputFormat{
	my $self=shift;
	$self->{'_outputFileType'}=shift // return $self->{'_outputFileType'};
}

sub inputFile{
	my $self=shift;
	$self->{'_inputFile'}=shift // return $self->{'_inputFile'};
}

sub conversionFile{
	my $self=shift;
	$self->{'_conversionFile'}=shift // return $self->{'_conversionFile'};
}

sub outputFile{
	my $self=shift;
	$self->{'_outputFile'}=shift // return $self->{'_outputFile'};
}

sub run{
	my($self)=shift;	
	
	$self->logger->info("Creating " . $self->outputFormat . " phylogeny file from " . $self->inputFile);

	my $inFH = IO::File->new('<' . $self->inputFile) // die "$!";

	#set the ordered list of names
	$self->nameOrderArray($self->_getHeaders($inFH->getline));

	#get the original data and transpose it for phylip format, send the number of columns we need
	$self->_printPhylipFormat(
		#uses Array::Transpose to return an arrayRef of the original data
		$self->_transposeArray(
			$self->_getOriginalData(
				$inFH, 
				scalar(@{$self->nameOrderArray})
			)
		)
	);

	$self->_printConversionInformation($self->nameOrderArray);
	$inFH->close();

	$self->logger->info("Finished");
}


=head2 _transposeArray

Takes the input array of arrays and transposes them.
Returns the transposed array of array ref.

=cut

sub _transposeArray{
	my $self = shift;
	my $arrayRef =shift;

	my @transposedArray;
	foreach my $row(@{$arrayRef}){
		for my $column(0..scalar(@{$row})){
			push(@{$transposedArray[$column]},$row->[$column]);
		}
	}
	return \@transposedArray;
}


=head2 _getOriginalData

Gathers the original data in an array of arrayRefs,
with each arrayRef being a row of data and returns a reference to the main array.
Only gathers the number of data columns needed for the generation of a phylogeny file.

=cut


sub _getOriginalData{
	my $self=shift;
	my $inFH = shift;
	my $numberOfColumns = shift;

	$self->logger->info("Number of data columns, including row ID: $numberOfColumns");

	my @originalData;
	while(my $line = $inFH->getline){
		my @la = split('\t',$line);
		my @neededColumns = splice(@la,1,($numberOfColumns -1) );

		push @originalData,\@neededColumns;
	}
	return \@originalData;
}


=head2 _getHeaders

Takes in the first row of data from the file, and creates a list of ordered headers.
This allows name conversion information to be generated as well as to identify the number
of data columns (ie. exclude the second and third set of columns with contig / location information)
that is stored in the panGenome and coreSnps files.
Returns an arrayRef to the ordered headers.

=cut

sub _getHeaders{
	my $self=shift;
	my $headerLine =shift;

	$headerLine =~ s/\R//g;
	my @la = split('\t',$headerLine);
	return \@la;
}


=head3 _printConversionInformation

Phylip format is limited to a 10-character name field.
In printing the Phylip format, we substitute numbers for names.
This creates a tab-delimited table that lists the conversion information.

=cut

sub _printConversionInformation{
	my $self=shift;
	my $arrayRef =shift;

	my $conversionFH = IO::File->new('>' . $self->conversionFile) or die "$!";
	$conversionFH->print(
		'#Name Conversion Information',
		"\n"
	);

	for my $i(1 .. (scalar(@{$arrayRef})-1)){
		$conversionFH->print(
			($i),
			"\t",
			$arrayRef->[$i],
			"\n"
		);
	}	
}

sub _printPhylipFormat{
	my($self)=shift;
	my $transposedArray = shift;

	my @names;
	my $numberOfGenomes = (scalar(@{$self->nameOrderArray})-1);

	my $outputFH=IO::File->new('>' . $self->outputFile) or die "$!";
	
	$outputFH->print(
		$numberOfGenomes,
		' ',
		scalar( @{$transposedArray->[0]} ),
		"\n" 
	);
	
	for my $i(0 .. ($numberOfGenomes -1)){
		$outputFH->print(
			($i+1),
			' 'x(10-length($i+1)),
			join('', @{$transposedArray->[$i]} ) . "\n"
		);
	}
	$outputFH->close();
}


1;




