#!/usr/bin/perl
package MSA::BlastBased::SNPFinder;

use FindBin;
use lib "$FindBin::Bin/../../";
use FileInteraction::Fasta::SequenceName;
use Log::Log4perl;

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_initialize(@_);
    return $self;
}

#class vars
sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

sub startHash{
	my $self=shift;
	$self->{'_startHash'} = shift // return $self->{'_startHash'};
}

sub sequenceHash{
	my $self=shift;
	$self->{'_sequenceHash'} = shift // return $self->{'_sequenceHash'};
}

sub allowableChars{
	my $self=shift;
	$self->{'_allowableChars'} = shift // return $self->{'_allowableChars'};
}

sub queryNameOrderHash{
	my $self=shift;
	$self->{'_queryNameOrderHash'}=shift // return $self->{'_queryNameOrderHash'};
}

sub firstSwitch{
	my $self=shift;
	$self->{'_firstSwitch'}=shift // return $self->{'_firstSwitch'};
}
#methods

sub _initialize{
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
	
	#logging
	$self->logger(Log::Log4perl->get_logger());
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
			
			if(defined $singleLine[0]){
				push @returnArray, (join("\t",@singleLine) . "\n");
			}
			else{
				$self->logger->debug("Single line empty");
			}
			
		}
		return \@returnArray;
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
		my $seqName = FileInteraction::Fasta::SequenceName->new($seq);
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
