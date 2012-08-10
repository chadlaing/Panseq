#!/usr/bin/perl

package TreeBuilding::PhylogenyFileCreator;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use IO::File;
use Log::Log4perl;
use parent 'FileInteraction::FlexiblePrinter';

sub new {
	my $class = shift;
	my $self= $class->SUPER::new(@_); #this calls FileInteraction::FlexiblePrinter::_initialize, so need to call the SUPER::_initialize
	return $self;
}

sub nameOrderArray{
	my $self=shift;
	$self->{'_nameOrderArray'}=shift // return $self->{'_nameOrderArray'};
}

sub tableHash{
	my $self=shift;
	$self->{'_tableHash'}=shift // return $self->{'_tableHash'};
}

sub isHeader{
	my $self=shift;
	$self->{'_isHeader'}=shift // return $self->{'_isHeader'};
}

sub phylipInfoFH{
	my $self=shift;
	$self->{'_phylipInfoFH'}=shift // return $self->{'_phylipInfoFH'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

#methods
sub _initialize{
	my $self=shift;
	
	#logging
	$self->logger(Log::Log4perl->get_logger());

	#inheritance
	$self->SUPER::_initialize(@_);

}

sub createFile{
	my($self)=shift;
	
	if(@_){
		my $type=shift;
		my $tabFile=shift;
		
		$self->logger->info("Creating $type phylogeny file from $tabFile");
		
		#deal with header/no-header files
		if(@_){
			$self->isHeader(1);
		}
		else{
			$self->isHeader(0);
		}
		
		$self->gatherTableInHash($tabFile) if (!defined $self->tableHash);
		
		if($type eq 'nexus'){
			$self->printNexusFormat();
		}
		elsif($type eq 'phylip'){
			$self->printPhylipFormat($tabFile);
		}
	}
	else{
		print STDERR "Nothing sent to createFile!\n";
		exit(1);
	}
}


sub printPhylipFormat{
	my($self)=shift;
	
	my @names;	
	
	$self->print(
		scalar(keys %{$self->tableHash}),
		' ',
		length($self->tableHash->{'1'}),
		"\n" 
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->tableHash});$i++){
		$self->print(
			$i,
			' 'x(10-length($i)),
			$self->tableHash->{$i} . "\n"
		);
		
		my $name;
		if(defined $self->nameOrderArray && defined $self->nameOrderArray->[$i-1]){
			$name = $self->nameOrderArray->[$i-1];
		}
		else{
			$name=$i;
		}
		push @names, $name;		
	}
	
	
	#print out the info conversion
	$self->outputFH($self->phylipInfoFH); #if phylipInfoFH is undefined, defaults to STDOUT (see FileInteraction::FlexiblePrinter.pm)
	$self->print(
		'#Name Conversion Information',
		"\n"
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->tableHash});$i++){
		$self->print(
			$i,
			"\t",
			$names[$i-1],
			"\n"
		);
	}	
}

sub printNexusFormat{
	my($self)=shift;
	
	my @names;
	
	$self->print(
		'#NEXUS' . "\n",
		'BEGIN Taxa;' . "\n",
		'DIMENSIONS ntax=' . scalar(keys %{$self->tableHash}) . ";\n",
		'TAXLABELS' . "\n"
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->tableHash});$i++){
		
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
		scalar(keys %{$self->tableHash}),
		' nchar=',
		length($self->tableHash->{1}) . ';' . "\n",
		'FORMAT datatype=dna symbols="ATGC" missing=? gap=-;' . "\n",
		'Matrix' . "\n",
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->tableHash});$i++){
		$self->print(
			$names[$i-1],
			"\t",
			$self->tableHash->{$i} . "\n"
		);
	}
	
	$self->print(
		';' . "\n" . 'End;'
	)
}

sub gatherTableInHash{
	my($self)=shift;
	
	if(@_){
		my $inFileName=shift;
		my $inFile = IO::File->new('<' . $inFileName) or die "$!";
		
		my $count=0;
		while(my $line = $inFile->getline){
			$count++;
			if($count==1 && $self->isHeader){
				$self->createNameOrderArray($line);
				next;
			}
			$line =~ s/\R//g;
			my $la = [split('\t',$line)];
			$self->addToTableHash($la);
		}		
		$inFile->close();
	}
	else{
		print "No file specified in gatherTableInHash!\n";
		exit(1);
	}
}

sub createNameOrderArray{
	my($self)=shift;
	
	if(@_){
		my $line=shift;
		$line =~ s/\R//g;
		my $la = [split('\t',$line)];
		
		my @orderedArray;
		for(my $i=1; $i<scalar(@{$la}); $i++){
			push @orderedArray, $la->[$i];
		}
		$self->nameOrderArray(\@orderedArray);
	}
	else{
		print STDERR "Nothing sent to createNameOrderArray!\n";
		exit(1);
	}
}


sub addToTableHash{
	my($self)=shift;
	
	if(@_){
		my $la=shift;
		
		 unless(defined $self->tableHash){
		 	my %tempHash;
		 	for(my $i=1; $i<scalar(@{$la}); $i++){
		 			$tempHash{$i}='';		
		 	}
		 	$self->tableHash(\%tempHash);
		 }
		
		for(my $i=1; $i<scalar(@{$la}); $i++){
			$self->tableHash->{$i} .= $la->[$i];			
		}
	}
	else{
		print STDERR "Nothing sent to addToTableHash!\n";
		exit(1);
	}
}

1;




