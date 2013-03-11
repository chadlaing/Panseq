#!/usr/bin/env perl

package Modules::Phylogeny::PhylogenyFileCreator;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use IO::File;
use Log::Log4perl;
use Role::Tiny::With;

with 'Roles::FlexiblePrinter';

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

    $self->logger->debug("Logger initialized in Modules::PanGenome::PanGenome");  

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
	$self->_tableHash({});
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

sub _tableHash{
	my $self=shift;
	$self->{'__tableHash'}=shift // return $self->{'__tableHash'};
}

sub outputFileType{
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

sub run{
	my($self)=shift;	
	
	$self->logger->info("Creating " . $self->outputFileType . " phylogeny file from " . $self->inputFile);

	$self->_tableHash($self->_gatherTableInHash);
	
	if($self->outputFileType eq 'nexus'){
		$self->printNexusFormat();
	}
	elsif($self->outputFileType eq 'phylip'){
		$self->printPhylipFormat();
	}

	if($self->conversionFile){
		$self->_printConversionInformation();
	}
}

=head3 _printConversionInformation

Phylip format is limited to a 10-character name field.
In printing the Phylip format, we substitute numbers for names.
This creates a tab-delimited table that lists the conversion information.

=cut

sub _printConversionInformation{
	my $self=shift;

	my $conversionFH = IO::File->new('>' . $self->conversionFile) or die "$!";
	$self->outputFH($conversionFH);
	$self->printOut(
		'#Name Conversion Information',
		"\n"
	);

	for my $i(1 .. (scalar(@{$self->nameOrderArray})-1)){
		$self->printOut(
			($i),
			"\t",
			$self->nameOrderArray->[$i],
			"\n"
		);
	}	
}

sub printPhylipFormat{
	my($self)=shift;
	
	my @names;
	my $numberOfGenomes = (scalar(@{$self->nameOrderArray})-1);
	
	$self->printOut(
		$numberOfGenomes,
		' ',
		length($self->_tableHash->{'1'}),
		"\n" 
	);
	
	for my $i(1 .. $numberOfGenomes){
		$self->printOut(
			($i),
			' 'x(10-length($i)),
			$self->_tableHash->{$i} . "\n"
		);
	}
}

sub printNexusFormat{
	my($self)=shift;
	
	my @names;
	
	$self->print(
		'#NEXUS' . "\n",
		'BEGIN Taxa;' . "\n",
		'DIMENSIONS ntax=' . scalar(keys %{$self->_tableHash}) . ";\n",
		'TAXLABELS' . "\n"
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->_tableHash});$i++){
		
		my $name;
		if(defined $self->nameOrderArray && defined $self->nameOrderArray->[$i-1]){
			$name = $self->nameOrderArray->[$i-1];
		}
		else{
			$name=$i;
		}
		push @names, $name;
		
		$self->print(
			$name,
			"\n"
		);
	}
	
	$self->print(
		'BEGIN data;' . "\n",
		'DIMENSIONS ntax=',
		scalar(keys %{$self->_tableHash}),
		' nchar=',
		length($self->_tableHash->{1}) . ';' . "\n",
		'FORMAT datatype=dna symbols="ATGC" missing=? gap=-;' . "\n",
		'Matrix' . "\n",
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->_tableHash});$i++){
		$self->print(
			$names[$i-1],
			"\t",
			$self->_tableHash->{$i} . "\n"
		);
	}
	
	$self->print(
		';' . "\n" . 'End;'
	)
}

sub _gatherTableInHash{
	my($self)=shift;

	my $inFile = IO::File->new('<' . $self->inputFile) or die "$!";
	
	my %tableHash;
	while(my $line = $inFile->getline){
		$line =~ s/\R//g;
		my @la = split('\t',$line);

		if($. == 1){
			$self->nameOrderArray(\@la);
			next;
		}

		for my $position(0..(scalar(@la)-1)){
			my $datum = $la[$position];
			$tableHash{$position} .=$datum;
		}
	}		
	$inFile->close();
	return \%tableHash;
}

1;




