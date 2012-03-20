#!/usr/bin/perl
package SNPFinder;

use FindBin::libs;
use FileInteraction::Fasta::SequenceName;

use Object::Tiny::RW qw{
	startHash
	sequenceHash
	allowableChars
	queryNameOrderHash
	firstSwitch
};

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->initialize(@_);
    return $self;
}

sub initialize{
	my($self)=shift;
	
	if(scalar(@_)==2){
		my $snpType=shift;
		my $orderHash=shift;
		
		$self->queryNameOrderHash($orderHash);
		
		if($snpType eq 'nucleotide'){
			$self->allowableChars('ACTGactg');
		}
		elsif($snpType eq 'protein'){
			$self->allowableChars('ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv');
		}
		else{
			print STDERR "Incorrect snpType set in SNPFinder!\n",
				"Allowable values are nucleotide and protein.\n";
			exit(1);
		}

		$self->firstSwitch('on');
	}
	else{
		print "SNPFinder requires a string of allowable characters to initialize!\n";
		exit(1);
	}
}

sub findSNPs{
	my($self)=shift;
	
	if(@_){
		#will return an array of SNP lines
		
		my $arrayRef=shift;
		my $countNumber=shift;
		
		#hash of startbp positions
		if(@_){
			my $startBpHashRef=shift;
			$self->startHash($startBpHashRef);
		}
		
		my $currentName;
		my %sequenceHash;
		my @arrayOfSNPHashRefs;
		my @returnArray;
		my @singleLine;
		
		#gather sequence
		foreach my $line(@{$arrayRef}){
			$line =~ s/[\n\f\r]//g;
			if($line =~ /^>(.+)/){
				$currentName=$1;
			}
			else{
				$sequenceHash{$currentName}.=$line;
			}
		}
		$self->sequenceHash(\%sequenceHash);
		
		#look through for SNPs
		my $seqLength=0;
		
		#get length
		foreach my $key(keys %sequenceHash){
			$seqLength = length($sequenceHash{$key});
			last;
		}
		
		for(my $i=0; $i< $seqLength; $i++){
			my $hashRef = $self->getHashRefOfSNP($i);
			
			if(defined $hashRef){
				push @arrayOfSNPHashRefs, $hashRef;
			}
		}	
		
		#create tab-delimited SNP line	
		foreach my $snpHashRef(@arrayOfSNPHashRefs){
			#get correct order
			my $infoLine;			
			foreach my $allQueryName(keys %{$self->queryNameOrderHash}){				
				if(defined $snpHashRef->{$allQueryName}){
					$singleLine[$self->queryNameOrderHash->{$allQueryName}]=$snpHashRef->{$allQueryName}->{'char'};
					$infoLine .= ($snpHashRef->{$allQueryName}->{'header'} . ',');
				}
				else{
					$singleLine[$self->queryNameOrderHash->{$allQueryName}]='-';
				}
			}
			$singleLine[0]=$infoLine;

			push @returnArray, (join("\t",@singleLine) . "\n");
		}
		return @returnArray;
	}
	else{
		print STDERR "No arrayRef of fasta sequences sent to initialize SNPFinder object!\n";
		exit(1);
	}
}

sub getHashRefOfSNP{
	my($self)=shift;
	
	my $position=shift;
	my $SNPHash={};
	my %charHash;
	my $allowableChars=$self->allowableChars;
	
	foreach my $seq(keys %{$self->sequenceHash}){
		my $seqName = SequenceName->new($seq);
		my $currChar = substr($self->sequenceHash->{$seq},$position,1);
		
		my $header; 
		if($self->firstSwitch eq 'on'){
			$header = '>' . (join('',@{$seqName->arrayOfHeaders}) . 'SNP@BP:');
		}
				
		if(defined $self->startHash){
			if(defined $self->startHash->{$seq}){
				$header .=($self->startHash->{$seq}+$position);
			}		
		}
		
		$SNPHash->{$seqName->name}->{'header'}= $header;
		$SNPHash->{$seqName->name}->{'char'}=$currChar;
		next unless $allowableChars =~ /$currChar/;
		$charHash{$currChar}=1;		
	}
	
	if(scalar(keys %charHash)>1){
		$self->firstSwitch('off');
		return $SNPHash;
	}
	else{
		return undef;
	}	
}

1;
