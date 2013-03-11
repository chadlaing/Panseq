#!/usr/bin/env perl
package Modules::PanGenome::TableExtractor;


use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Role::Tiny;
use Log::Log4perl;

with 'Roles::FlexiblePrinter';


=head1 Description

Takes the output tables of Panseq (core_snps.txt) and (pan_genome.txt) and produces new tables
based on the specified parameters.
For example, take only SNP lines that contain fewer than X '-', and/or that have
at least 2 SNP characters present at least Y times.
Or, for pan_genome data, produce a binary table based on a new %ID cutoff value.

=cut

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

=head2 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}


=head3 _initialize

Initializes the logger.
Assigns all values to class variables.
Anything else that the _initialize function does.

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 
    $self->logger->debug("Logger initialized in Modules::PanGenome::TableExtractor");  

    my %params = @_;

    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::PanGenome::TableExtractor");
		}
	}	
}

sub inputFile{
	my $self=shift;
	$self->{'_inputFile'}=shift // return $self->{'_inputFile'};
}

sub percentId{
	my $self=shift;
	$self->{'_percentId'}=shift // return $self->{'_percentId'};
}

sub minimumCharacters{
	my $self=shift;
	$self->{'_minimumDifferentCharacters'}=shift // return $self->{'_minimumDifferentCharacters'};
}

sub maximumMissing{
	my $self=shift;
	$self->{'_maximumMissing'}=shift // return $self->{'_maximumMissing'};
}

sub minimumPresent{
	my $self=shift;
	$self->{'_minimumPresent'}=shift // return $self->{'_minimumPresent'};
}

sub tableType{
	my $self=shift;
	$self->{'_tableType'}=shift // return $self->{'_tableType'};
}

=head2 run

Runs the program. Calls _extractCoreTable or _extractPanTable
depending on the tableType value.

=cut

sub run{
	my $self =shift;
	$self->logger->info("Start");

	unless(defined $self->tableType){
		$self->logger->fatal("tableType required in TableExtractor");
	}

	my $inFH = IO::File->new('<' . $self->inputFile) or die "Could not open file $!";

	if($self->tableType eq 'core'){
		$self->_extractCoreTable($inFH);
	}
	elsif($self->tableType eq 'pan'){
		$self->_extractPanTable($inFH);
	}
	else{
		$self->logger->fatal("\nIncorrect tableType of " . $self->tableType
			. ".\n The values \'core\' and \'pan\' are acceptable."
		);
		exit(1);
	}
	$self->logger->info("End");
}


=head2 _extractCoreTable

Takes a tab-delimited file of SNPs, and prints only the lines that contain fewer than the
user specified number of gaps, and at least X of each SNP character.

=cut

sub _extractCoreTable{
	my $self=shift;
	my $inFH = shift;

	if((defined $self->minimumPresent && defined $self->maximumMissing) || 
		(!defined $self->minimumPresent && !defined $self->maximumMissing)){
			$self->logger->fatal("\nExactly one of minimumPresent or maximumMissing must be defined at one time");
			exit(1);
	}

	my $numberOfColums=0;
	my $stringToMatch;
	while(my $line = $inFH->getline){
		if($inFH->input_line_number == 1){
			$numberOfColums = $self->_determineNumberOfColumns($line);
			$self->logger->info("Number of columns: $numberOfColums");

			unless(defined $self->minimumPresent){
				$self->minimumPresent($numberOfColums - $self->maximumMissing);
				$self->logger->info("minimumPresent: " . $self->minimumPresent);
			}

			#determine the number of \t\w instances to match
			$stringToMatch = $self->_createStringToMatch($numberOfColums);
			next;
		}

		my $dataLine;		
		if($line =~ m/($stringToMatch)/){
			$dataLine = $1;
		}
		else{
			$self->logger->info("Dataline not defined!");
			exit;
		}	

		my $numberPresent = () = $dataLine =~ m/\t[^-]/g;	

		if($numberPresent >= $self->minimumPresent){
			#use Roles::FlexiblePrinter->printOut (default STDOUT)
			$self->printOut($line . "\n");
		}
	}
	$inFH->close();
}

=head2 _createStringToMatch

Based on the number of columns, create a perl string to use for a regex match

=cut

sub _createStringToMatch{
	my $self = shift;
	my $numberOfColumns = shift;

	my $matchString='^\S+';
	for(1..$numberOfColumns){
		$matchString .= '\t[\w\-]';
	}
	return qr/$matchString/;
}

=head2 _determineNumberOfColumns

Given the output table first line, return a count of how many data columns

=cut

sub _determineNumberOfColumns{
	my $self=shift;
	my $line =shift;

	my $count = () = $line =~ m/\t\S+/g;
	return $count;
}


=head2 _extractPanTable

Given the input file of tab-delimited percentIds, and a percentId cutoff, prints out the lines as binary,
specifying 1 as greater than or equal the percentId cutoff and 0 as below.

=cut

sub _extractPanTable{
	my $self =shift;
	my $inFH = shift;

	unless(defined $self->percentId){
		$self->logger->fatal("\nperdentId is required in a 'pan' run\n");
		exit(1);
	}

	my $numberOfColums;
	while(my $line = $inFH->getline){
		if($inFH->input_line_number == 1){
			$numberOfColums = $self->_determineNumberOfColumns($line);
			$self->logger->info("Number of columns: $numberOfColums");
			next;
		}

		$self->printOut($self->_getBinaryLine($line,$numberOfColums));
	}
}


=head2 _getBinaryLine

Based on the percentId cutoff-value, given a line of tab-delimited percent-Ids, returns a binary 
representation of the line of percentIds present in the input line.

=cut


sub _getBinaryLine{
	my $self=shift;
	my $line =shift;
	my $numberOfColums=shift;

	my $count=0;
	while($line =~ m/\t(\d)/gc){
		my $percentId = $1;

		if($count > $numberOfColums){
			last;
		}	

		if($percentId >= $self->percentId){
			$line =~ s/$percentId/1/;
		}
		else{
			$line =~ s/$percentId/0/;
		}
		$count++;
	}
	return $line;	
}

1;



