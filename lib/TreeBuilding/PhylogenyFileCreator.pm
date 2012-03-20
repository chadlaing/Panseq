#!/usr/bin/perl

package PhylogenyFileCreator;

use FindBin::libs;
use IO::File;
use FileInteraction::LinePackage;
use FileInteraction::FlexiblePrinter;

our @ISA = qw{FlexiblePrinter}; #includes outputFilehandle and printOut

#object creation
use Object::Tiny::RW qw{
	nameOrderArray
	tableHash
	isHeader
	phylipInfoFH
};


sub createFile{
	my($self)=shift;
	
	if(@_){
		my $type=shift;
		my $tabFile=shift;
		
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
	
	$self->printOut(
		scalar(keys %{$self->tableHash}),
		' ',
		length($self->tableHash->{'1'}),
		"\n" 
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->tableHash});$i++){
		$self->printOut(
			$i,
			' 'x(10-length($i)),
			$self->tableHash->{$i} . "\n"
		);
		
		if(defined $self->nameOrderArray && defined $self->nameOrderArray->[$i-1]){
			$name = $self->nameOrderArray->[$i-1];
		}
		else{
			$name=$i;
		}
		push @names, $name;		
	}
	
	
	#print out the info conversion
	$self->outputFilehandle($self->phylipInfoFH); #if phylipInfoFH is undefined, defualts to STDOUT (see FlexiblePrinter.pm)
	$self->printOut(
		'#Name Conversion Information',
		"\n"
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->tableHash});$i++){
		$self->printOut(
			$i,
			":\t",
			$names[$i-1],
			"\n"
		);
	}	
}

sub printNexusFormat{
	my($self)=shift;
	
	my @names;
	
	$self->printOut(
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
		
		$self->printOut(
			$name,
			"\n"
		);
	}
	
	$self->printOut(
		'BEGIN data;' . "\n",
		'DIMENSIONS ntax=',
		scalar(keys %{$self->tableHash}),
		' nchar=',
		length($self->tableHash->{1}) . ';' . "\n",
		'FORMAT datatype=dna symbols="ATGC" missing=? gap=-;' . "\n",
		'Matrix' . "\n",
	);
	
	for(my $i=1; $i<=scalar(keys %{$self->tableHash});$i++){
		$self->printOut(
			$names[$i-1],
			"\t",
			$self->tableHash->{$i} . "\n"
		);
	}
	
	$self->printOut(
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
			my $la = LinePackage->new($line);
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
		my $la = LinePackage->new($line);
		
		my @orderedArray;
		for(my $i=1; $i<scalar(@{$la->lineArray}); $i++){
			push @orderedArray, $la->lineArray->[$i];
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
			for(my $i=1; $i<scalar(@{$la->lineArray}); $i++){
					$tempHash{$i}='';		
			}
			$self->tableHash(\%tempHash);
		}
		
		for(my $i=1; $i<scalar(@{$la->lineArray}); $i++){
			$self->tableHash->{$i} .= $la->lineArray->[$i];			
		}
	}
	else{
		print STDERR "Nothing sent to addToTableHash!\n";
		exit(1);
	}
}

1;




