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
	if($self->tableType eq 'core'){
		$self->_extractCoreTable();
	}
	elsif($self->tableType eq 'pan'){
		$self->_extractPanTable();
	}
	else{
		$self->logger->fatal("\nIncorrect tableType of " . $self->tableType
			. ".\n The values \'core\' and \'pan\' are acceptable."
		);
		exit(1);
	}
	$self->logger->info("End");
}


sub _extractCoreTable{
	my $self=shift;

	if((defined $self->minimumPresent && defined $self->maximumMissing) || 
		(!defined $self->minimumPresent && !defined $self->maximumMissing)){
			$self->logger->fatal("\nExactly one of minimumPresent or maximumMissing must be defined at one time");
			exit(1);
	}

	my $inFH = IO::File->new('<' . $self->inputFile) or die "Could not open file $!";

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
			#$self->logger->info("stringToMatch: $stringToMatch");
			next;
		}

		my $dataLine;		
		if($line =~ m/($stringToMatch)/){
			$dataLine = $1;
			#$self->logger->info("dataline: $dataLine");
		}
		else{
			$self->logger->info("Dataline not defined!");
			exit;
		}	

		my $numberPresent = () = $dataLine =~ m/\t[^-]/g;	
		#$self->logger->info("Number present: $numberPresent");

		if($numberPresent >= $self->minimumPresent){
			#use Roles::FlexiblePrinter->printOut (default STDOUT)
			#$self->printOut($dataLine . "\n");
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

1;



