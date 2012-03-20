#!/usr/bin/perl

#written by Chad Laing; last updated October 05, 2011

#Pan-genome sequence analysis using Panseq: an online tool for the rapid analysis of core and accessory genomic regions
#Chad Laing, Cody Buchanan, Eduardo Taboada, Yongxiang Zhang, Andrew Kropinski, Andre Villegas, James E Thomas and Victor PJ Gannon
#BMC Bioinformatics 2010, 11:461
#http://www.biomedcentral.com/1471-2105/11/461

package LociSelector;

#usage
# my $obj = $LociSelector->new();
# $obj->getBestLoci(<file name (mandatory)>,<loci number or 'best (mandatory)>,<output filehandle (optional, STDOUT defualt)>);

use FindBin::libs;
use IO::File;
use FileInteraction::LinePackage;
use FileInteraction::FlexiblePrinter;

our @ISA = qw{FlexiblePrinter}; #includes outputFilehandle and printOut

#object creation
use Object::Tiny::RW qw{
	lociNumber
	PODHash
	characterHash
	fingerprintHash
	missingCharsHash
	orderedLociHash
	outputLociHash
	toContinue
	currentFingerprints
	strainNames
	maskedPODLociHash
};

sub new{
	my($class)  = shift;

    my $self= {};
    bless ($self, $class);
    
    #defaults the characters not included in POD or fingerprint calculations to '? and '-'
    #by calling ->new(<list of chars>), one can personalize the missing charcters
    if(@_){
    	$self->setMissingChars(@_);
    }
    else{
    	$self->setMissingChars('-','?');
    }
    
    return $self;
}

#methods
sub setMissingChars{
	my($self)=shift;
	
	if(@_){
		foreach(@_){
			$self->addToMissingCharsHash($_);
		}
	}
	else{
		print STDERR "nothing sent to setMissingChars\n";
		exit(1);
	}
}

sub addToMissingCharsHash{
	my($self)=shift;
	
	#takes arguments as key
	
	if(scalar(@_)==1){
		if(defined $self->missingCharsHash){
			$self->missingCharsHash->{$_[0]}=1;
		}
		else{
			my %tempHash;
			$tempHash{$_[0]}=1;
			$self->missingCharsHash(\%tempHash);
		}
	}
	else{
		print STDERR "wrong number of arguments to addToMissingCharsHash\n";
		exit(1);
	}
}

sub getBestLoci{
	my($self)=shift;
	
	if(scalar(@_) >=2){
		my $inFile = IO::File->new('<'. $_[0]) or die "$!";
		$self->lociNumber($_[1]);
		
		if($_[2]){
			$self->outputFilehandle($_[2]);
		}
		
		$self->getAllLoci($inFile);
		$self->runLociSelector();
		$inFile->close();
	}
	else{
		print STDERR "incorrect number of arguments to getBestLoci!\n";
		exit(1);
	}
}

sub runLociSelector{
	my($self)=shift;
	
	#takes no arguments
	#algorithm:
	#choose best seed locus based on #fingerprints, then #POD
	#add subsequent loci until:
	#	a)all unique fingerprints have been exhausted or
	#	b)the # of fingerprints specified by the user has been reached
	#create fingerprints for each strain
	
	my $toContinue=1;
	my $exhaustedFingerprints=0;
	while($toContinue==1){
		my @locusInfo; #stores [0]=locus, [1]=POD
		
		if($exhaustedFingerprints){
			$locusInfo[0]=0;
		}
		else{
			@locusInfo=$self->chooseNextBestFingerprintLocus();
			$self->updateFingerprintHash($locusInfo[0]) if $locusInfo[0] ne '0';
			#if $locus eq 0, then unique fingerprints have been exhausted
		}
		
		#throw switches
		if($locusInfo[0] eq '0'){
			$exhaustedFingerprints=1;
			$toContinue = 0 if $self->lociNumber eq 'best';
		}
		
		#get next best POD if unique fingerprints exhausted
		if($toContinue && ($locusInfo[0] eq '0')){
			@locusInfo = $self->chooseNextBestPODLocus(); #this returns [0]=locus, [1]=POD	
			$toContinue =0 if $locusInfo[0] eq '0';
			$self->maskLociPair($locusInfo[0]);		
		}
		#update the output
		$self->addToOutputLociHash(@locusInfo) if $locusInfo[0] ne '0';
		
		#stop once the required loci have been gathered
		$toContinue = 0 if($self->lociNumber ne 'best' && defined $self->outputLociHash && scalar(keys %{$self->outputLociHash})==$self->lociNumber);
		
		#check to ensure that after the current locus is removed, there are still loci available
		$toContinue = 0 if(scalar(keys %{$self->characterHash}) ==1);
		
		#remove loci from future possible choices
		if($toContinue){		
			$self->removeLocusFromCharacterHash($locusInfo[0]);			
		}		
	}
	$self->printResults();
}

sub printResults{
	my($self)=shift;
	
	$self->printOut(localtime() . "\n",
	'FINAL VALUES' . "\n" . "==========" . "\n",
	'Target Loci Number:' . $self->lociNumber . "\n");
	$self->printOut('No. Of Loci Selected: ' . scalar(keys %{$self->outputLociHash}) . "\n") if $self->lociNumber eq 'best';
	$self->printOut('Fingerprints: ' . $self->currentFingerprints . "\n");
	
	my %outputHash;
	foreach my $locus(keys %{$self->outputLociHash}){
		my $outputLine = $locus . "\t" .
			 'POD: ' . $self->outputLociHash->{$locus}->{'POD'} . "\t";
		
		for(my $i=1; $i <  scalar (@{$self->outputLociHash->{$locus}->{'characters'}->lineArray}); $i++){
			$outputLine .= $self->outputLociHash->{$locus}->{'characters'}->lineArray->[$i];
		}
		$outputHash{$self->outputLociHash->{$locus}->{'order'}}=$outputLine;
	}
	
	#print out ordered output
	for(my $i=1; $i <= scalar(keys %outputHash); $i++){
		$self->printOut($outputHash{$i} . "\n");
	}	
}

sub maskLociPair{
	my($self)=shift;
	
	if(@_){
		my $locus=shift;
		
		foreach my $pair(%{$self->PODHash}){
			foreach my $pairLocus(%{$self->PODHash->{$pair}}){
				if($locus eq $pairLocus){
					$self->addToMaskedPODLociHash($pair);
					last;
				}
			}
		}
	}
	else{
		print STDERR "nothing sent to maskLocus\n";
		exit(1);
	}	
}

sub chooseNextBestPODLocus{
	my($self)=shift;
	
	#if send an arrayRef of loci, only choose those, otherwise, use all loci
	my $lociToTestRef;
	if(@_){
		$lociToTestRef=shift;
	}
	else{
		my @remainingLoci = keys %{$self->characterHash};
		$lociToTestRef = \@remainingLoci;
	}
	
	my @bestLocusWithPOD = $self->getHighestNonMaskedPODLocus($lociToTestRef);

	#if everything is masked, reset and try again
	if($bestLocusWithPOD[0] eq '0'){
		$self->resetMaskedPODLoci();
		@bestLocusWithPOD = $self->getHighestNonMaskedPODLocus($lociToTestRef);
	}
	else{
		return @bestLocusWithPOD;
	}
	
	
}

sub getHighestNonMaskedPODLocus{
	my($self)=shift;
	
	if(@_){
		my $arrayRef=shift;		
		my $highestPOD=0;
		my $highestLocus=0;
		
		#get top POD locus	
		my $PODRef = $self->calculateCurrentPOD($arrayRef);
		foreach my $locus(keys %{$PODRef}){
			if($PODRef->{$locus} > $highestPOD){
				$highestPOD = $PODRef->{$locus};
				$highestLocus = $locus;
			}
		}		
		return($highestLocus, $highestPOD);
	}
	else{
		print STDERR "nothing sent to getHighestNonMaskedPODLocus\n";
		exit(1);
	}
}

sub calculateCurrentPOD{
	my($self)=shift;
	
	#sent an arrayRef of loci we would like POD for, given the current masking parameters	
	#returns a hashRef for $hash{locus}=POD
	if(scalar(@_)==1){
		my $arrayRef=shift;
		my %PODHash;

		#create a hash of loci we are interested in
		my %lociWeCareAbout;
		foreach(@{$arrayRef}){
			$lociWeCareAbout{$_}=1;
		}
		
		
		#generate a count of POD for each allowable locus we care about
		foreach my $locusPair(keys %{$self->PODHash}){
			next if((defined $self->maskedPODLociHash) && (defined $self->maskedPODLociHash->{$locusPair}));
			
			foreach my $locus(keys %{$self->PODHash->{$locusPair}}){
				next unless defined $lociWeCareAbout{$locus};
			
				if(defined $PODHash{$locus}){
					$PODHash{$locus}++;
				}
				else{
					$PODHash{$locus}=1;
				}
			}
		}
		
		#if there are no allowable hits for any of the loci, return 0 for POD
		unless (%PODHash){
			foreach(keys %lociWeCareAbout){
				$PODHash{$_}=0;
			}
		}
	
		return \%PODHash;
	}
	else{
		print STDERR "nothing sent to calculateCurrentPOD\n";
		exit(1);
	}
}


sub removeLocusFromCharacterHash{
	my($self)=shift;
	
	if(@_){
		my $locus=shift;
		delete $self->characterHash->{$locus};
	}
	else{
		print STDERR "nothing sent to removeLocusFromCharacterHash\n";
		exit(1);
	}
}

sub addToOutputLociHash{
	my($self)=shift;
	
	#takes in an array, not a ref, of [0]=locus, [1]=POD
	#automatically adds the lineArray object for the locus into the hash for output as {'characters'}
	#POD in hash is {'POD'} for each locus
	
	if(@_){
		my $newLocus=shift;
		my $POD=shift;
		
		if(!defined $self->outputLociHash){
			my %tempHash;
			$tempHash{$newLocus}->{'characters'}=$self->characterHash->{$newLocus};
			$tempHash{$newLocus}->{'POD'}=$POD;
			$tempHash{$newLocus}->{'order'}=1;
			$self->outputLociHash(\%tempHash);
		}
		else{
			die "cannot find locus: $newLocus!\n" if !defined $self->characterHash->{$newLocus};
			die "cannot find POD: $POD!\n" if !defined $POD;
			$self->outputLociHash->{$newLocus}->{'characters'}=$self->characterHash->{$newLocus};
			$self->outputLociHash->{$newLocus}->{'POD'}=$POD;
			$self->outputLociHash->{$newLocus}->{'order'}= scalar(keys %{$self->outputLociHash});
		}
	}
	else{
		print STDERR "nothing sent to addToOutputLociHash\n";
		exit(1);
	}
}

sub updateFingerprintHash{
	my($self)=shift;
	
	if(@_){
		my $locus=shift;
		
		for(my $i=1; $i< scalar(@{$self->characterHash->{$locus}->lineArray});$i++){
			my $currentFingerprint;

			if(defined $self->fingerprintHash && defined $self->fingerprintHash->{$i}){
				$currentFingerprint=$self->fingerprintHash->{$i} . $self->characterHash->{$locus}->lineArray->[$i];
			}
			else{
				$currentFingerprint=$self->characterHash->{$locus}->lineArray->[$i];
			}	
			#make sure that any missing chars are converted to '.' in the string
			foreach my $missingChar(keys %{$self->missingCharsHash}){
				$currentFingerprint =~ s/\Q$missingChar/./;
			}
					
			$self->addToFingerprintHash($i,$currentFingerprint);
		}		
	}
	else{
		print STDERR "no locus sent for updating!\n";
		exit(1);
	}
}

sub chooseNextBestFingerprintLocus{
	my($self)=shift;
	
	my $bestLociForFingerprintRef = $self->getBestLociForFingerprint;
	
	my @locusToReturnInfo;
	
	if((scalar @{$bestLociForFingerprintRef}) > 1){
		@locusToReturnInfo = $self->chooseNextBestPODLocus($bestLociForFingerprintRef);
	}
	elsif((scalar @{$bestLociForFingerprintRef}) == 1){
		$locusToReturnInfo[0] = $bestLociForFingerprintRef->[0];
		my $PODHashRef = $self->calculateCurrentPOD($bestLociForFingerprintRef);
		$locusToReturnInfo[1]=$PODHashRef->{$locusToReturnInfo[0]};
	}
	else{
		$locusToReturnInfo[0]=0;
	}	
	return @locusToReturnInfo;
}

sub resetMaskedPODLoci{
	my($self)=shift;
	my %emptyHash;
	$self->maskedPODLociHash(\%emptyHash);
}

sub getBestLociForFingerprint{
	my($self)=shift;
	
	if(defined $self->characterHash){
		my $currentFingerprints= ($self->currentFingerprints || 0);
		my @bestFingerprintLoci=(0);
		my $doesItMakeThingsBetter=0; #switch so that we only return loci which increase the current # of fingerprints
		
		foreach my $locus(keys %{$self->characterHash}){
			my $returnedNumFingerprints = $self->calculateFingerprint($locus);
			
			if($returnedNumFingerprints > $currentFingerprints){
				$doesItMakeThingsBetter=1;
				@bestFingerprintLoci=();
				push @bestFingerprintLoci, $locus;
				$currentFingerprints = $returnedNumFingerprints;
				$self->currentFingerprints($currentFingerprints);
			}
			elsif(($returnedNumFingerprints == $currentFingerprints) && ($doesItMakeThingsBetter==1)){
				push @bestFingerprintLoci, $locus;
			}
		}
		return \@bestFingerprintLoci;
	}
	else{
		@bestFingerprintLoci=(0);
		return \@bestFingerprintLoci;
	}
}

sub calculateFingerprint{
	my($self)=shift;
	
	if(@_){
		my $newLocus=shift;
		my @fingerArray;
		my $numberOfUniqueFingerprints=0;
		my @observedPrints;
	
		#build an array of fingerprints based on the selected loci
		#initialize with current fingerprint values
		if(defined $self->fingerprintHash){
			foreach(keys %{$self->fingerprintHash}){
				$fingerArray[$_]=$self->fingerprintHash->{$_};
			}
		}		

		for(my $i=1; $i< scalar(@{$self->characterHash->{$newLocus}->lineArray});$i++){
			my $char = $self->characterHash->{$newLocus}->lineArray->[$i];
			my @tempArray=($char);
			$char = '.' unless $self->isValidChars(\@tempArray);
			if(defined $fingerArray[$i]){
				$fingerArray[$i].=$char;
			}
			else{
				$fingerArray[$i]=$char;
			} 
		}	
		
		#count unique fingerprints
		ALLPRINTS: for(my $i=1; $i< scalar(@fingerArray); $i++){
			my $aPrint = $fingerArray[$i];
			unless(@observedPrints){
				push @observedPrints, $aPrint;
				$numberOfUniqueFingerprints++;
				next ALLPRINTS;
			}
			
			OBSERVEDPRINTS: foreach my $oPrint(@observedPrints){
				if(($aPrint =~ m/$oPrint/) || ($oPrint =~ m/$aPrint/)){
					next ALLPRINTS;
				}
				
			}
			$numberOfUniqueFingerprints++;
			push @observedPrints, $aPrint;			
		}
		
		return $numberOfUniqueFingerprints;
	}
	else{
		print STDERR "nothing sent to calculateFingerprints\n";
		exit(1);
	}
}

sub getAllLoci{
	my($self)=shift;
	
	if(@_){
		my $inFile=shift;
		
		while(my $line = $inFile->getline){
			my $la = LinePackage->new($line);
			my $locusName = $la->lineArray->[0];
			
			if($.==1){
				#initialize the fingerprint hash on first line
				$self->strainNames($la);
			}
			else{
				#add the character values for the locus to the characterHash
				$self->addToCharacterHash($locusName,$la);
				
				#calculate the POD
				my $allLociRef = $self->calculatePOD($la);				
				foreach my $differentPair(@{$allLociRef}){
					$self->addToPODHash($differentPair,$locusName);
				}				
			}
		}
	}
	else{
		print STDERR "nothing sent to getPODofAllLoci\n";
		exit(1);
	}
}

sub addToCharacterHash{
	my($self)=shift;
	
	#takes arguments as key,value
	
	if(scalar(@_)==2){
		if(defined $self->characterHash){
			$self->characterHash->{$_[0]}=$_[1];
		}
		else{
			my %tempHash;
			$tempHash{$_[0]}=$_[1];
			$self->characterHash(\%tempHash);
		}
	}
	else{
		print STDERR "wrong number of arguments to addToCharacterHash\n";
		exit(1);
	}
}

sub addToMaskedPODLociHash{
	my($self)=shift;
	
	#takes argument as the locus pair ($i$j) to mask
	
	if(scalar(@_)==1){
		my $maskedPair=shift;
		if(defined $self->maskedPODLociHash){
			$self->maskedPODLociHash->{$maskedPair}=1;
		}
		else{
			my %tempHash;
			$tempHash{$maskedPair}=1;
			$self->maskedPODLociHash(\%tempHash);
		}
	}
	else{
		print STDERR "wrong number of arguments to addToMaskedPODLociHash\n";
		exit(1);
	}
}

sub calculatePOD{
	my($self)=shift;
	
	if(@_){
		my $la=shift;		
		my @differencesArray;
		
		for(my $i=1; ($i+1)< scalar(@{$la->lineArray});$i++){			
			my $firstChar=$la->lineArray->[$i];
			
			for(my $j=($i+1); $j < scalar(@{$la->lineArray});$j++){
				my $secondChar=$la->lineArray->[$j];				
				if($self->isDifferentChars($firstChar,$secondChar)){
					push @differencesArray, "${i}vs${j}";
				}
			}		
		}
		return \@differencesArray;
	}
	else{
		print STDERR "nothing sent to calculatePOD\n";
		exit(1);
	}
}

sub addToPODHash{
	my($self)=shift;
	
	#takes arguments as key,value
	
	if(scalar(@_)==2){
		if(defined $self->PODHash){
			$self->PODHash->{$_[0]}->{$_[1]}=1;
		}
		else{
			my %tempHash;
			$tempHash{$_[0]}->{$_[1]}=1;
			$self->PODHash(\%tempHash);
		}
	}
	else{
		print STDERR "wrong number of arguments to addToPODHash\n";
		exit(1);
	}
}


sub isDifferentChars{
	my($self)=shift;
	
	if(scalar(@_)==2){
		
		if($self->isValidChars(\@_) && ($_[0] ne $_[1])){
			return 1;
		}
		else{
			return 0;
		}
	}
	else{
		print STDERR "incorrect number of arguments sent to isDifferentChars\n";
		exit(1);
	}
}

sub isValidChars{
	my($self)=shift;
	
	if(@_){
		my $arrayRef=shift;
		foreach(@{$arrayRef}){
			return 0 if defined $self->missingCharsHash->{$_};
		}
		return 1;
	}
	else{
		print STDERR "nothing sent to isValidChars\n";
		exit(1);
	}
}

sub addToFingerprintHash{
	my($self)=shift;
	
	#takes arguments as key,value
	
	if(scalar(@_)==2){
		if(defined $self->fingerprintHash){
			$self->fingerprintHash->{$_[0]}=$_[1];
		}
		else{
			my %tempHash;
			$tempHash{$_[0]}=$_[1];
			$self->fingerprintHash(\%tempHash);
		}
	}
	else{
		print STDERR "wrong number of arguments to addToFingerprintHash\n";
		exit(1);
	}
}

1;
