#!/usr/bin/perl

package CoreAccessoryProcessor;

use strict;
use warnings;
use diagnostics;
use Carp;
use FindBin::libs;
use MSA::BlastBased::BlastResultFactory;
use MSA::BlastBased::BlastResultObject;
use MSA::BlastBased::BlastHitObject;
use MSA::BlastBased::SNPFinder;
use Muscle::MuscleCmd;
use MSA::BlastBased::CoreAccessory;
use File::Temp;

our @ISA= qw/CoreAccessory/;

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_coreAccessoryProcessorInitialize(@_);
    return $self;
}

#class variables
sub _queryNameOrderHash{
	my $self=shift;
	$self->{'_CoreAccessoryProcessor_queryNameOrderHash'}=shift // return $self->{'_CoreAccessoryProcessor_queryNameOrderHash'};
}

sub _accessoryTempFile{
	my $self=shift;
	$self->{'_CoreAccessoryProcessor_accessoryTempFile'}=shift // return $self->{'_CoreAccessoryProcessor_accessoryTempFile'};
}


sub _coreTempFile{
	my $self=shift;
	$self->{'_CoreAccessoryProcessor_coreTempFile'}=shift // return $self->{'_CoreAccessoryProcessor_coreTempFile'};
}

#methods
sub _coreAccessoryProcessorInitialize{
	my $self=shift;
	
	#set values
	my $paramsRef=shift;
	$self->_baseDirectory($paramsRef->{'baseDirectory'}) // confess ('baseDirectory required in coreAccessoryProcessor');
	$self->_percentIdentityCutoff($paramsRef->{'percentIdentityCutoff'}) // confess ('percentIdentityCutoff required in coreAccessoryProcessor');
	$self->queryNameObjectHash($paramsRef->{'queryNameObjectHash'}) // confess ('queryNameObjectHash required in coreAccessoryProcessor');
	$self->_coreGenomeThreshold($paramsRef->{'coreGenomeThreshold'}) // confess ('coreGenomeThreshold required in coreAccessoryProcessor');
	$self->_accessoryType($paramsRef->{'accessoryType'}) // confess ('accessoryType required in coreAccessoryProcessor');
	$self->_snpType($paramsRef->{'snpType'}) // confess ('snpType required in coreAccessoryProcessor');		
	
	#inheritance
	$self->_coreAccessoryInitialize(@_);
	
	#set up queryNameOrderHash
	$self->_getQueryNameOrder(); 
}


sub _getQueryNameOrder{
	my($self)=shift;
	
	my %returnHash;
	my $order=1;
	
	$self->logger->debug("DEBUG:\tRef type of queryNameObjectHash" . ref($self->queryNameObjectHash));
	$self->logger->debug("DEBUG:\tScalar of queryNameObjectHash" . scalar keys %{$self->queryNameObjectHash});
	
	
	foreach my $name(keys %{$self->queryNameObjectHash}){
		$returnHash{$name}=$order;
		$order++;
	}
	
	$self->_queryNameOrderHash(\%returnHash);
}


sub _combineTempFiles{
	my($self)=shift;
	
	if(scalar(@_)==2){
		my $textToMatch=shift;
		my $outputFileName=shift;
		
		my @filesToCombine;
		my $fm = FileManipulation->new();
		my $namesRef = $fm->getFileNamesFromDirectory($self->_baseDirectory);
		
		foreach my $name(@{$namesRef}){
			if($name =~ /$textToMatch\d+/){
				push @filesToCombine, ($name);
			}
		}		
		
		if(scalar(@filesToCombine) > 0){			
			my $oFile = IO::File->new('>' . $outputFileName) or die "$!";
			
			#print out header names on combined files (tab-delimited names)
			my @outputArray;
			foreach my $name(keys %{$self->_queryNameOrderHash}){
				$outputArray[$self->_queryNameOrderHash->{$name}]=$name;
			}
			
			for(my $i=1; $i<scalar(@outputArray);$i++){
				print $oFile "\t" . $outputArray[$i];
			}
			print $oFile "\n";
			
			$fm->outputFilehandle($oFile);
			$fm->combineFilesIntoSingleFile(\@filesToCombine);		
			$oFile->close();	
			
			#remove temp files	
			foreach(@filesToCombine){
				#unlink;
			}
		}
	}
	else{
		print STDERR "Wrong number of arguments to combineTempFiles!\n";
		exit(1);
	}
}


sub processBlastXML{
	my($self)=shift;
	
	if(scalar(@_)==2){
		my $blastXMLFilesArrayRef=shift;
		my $numberOfCores=shift;
				
		my $forker = Parallel::ForkManager->new($numberOfCores);
		my $currentResult=0;
		
		$self->logger->info("INFO:\tProcessing " . scalar(@{$blastXMLFilesArrayRef}) . " Blast XML files");
		
		foreach my $xmlFile(@{$blastXMLFilesArrayRef}){
			$currentResult++;
			
			$self->logger->debug("DEBUG:\tProcessing xmlFile: $xmlFile");
						
			$forker->start and next;
				
				#create filehandle
				my $xmlFH = IO::File->new('<' . $xmlFile) or die "$!";
				my $xmler = BlastResultFactory->new($xmlFH);
								
				my @resultArray;
				while(my $result = $xmler->nextResult){
					push @resultArray, $result;
				}
				$self->_sendToProcessQueue(\@resultArray, $currentResult) if (scalar(@resultArray) >0);	
				$xmlFH->close();
			$forker->finish;
		}
		$forker->wait_all_children();
		
		$self->_combineResultTempFiles($numberOfCores);
		
		#process the combined files
	
		
	}#end if
	else{
		print STDERR "please specify a BLAST xml output file for processing!\n";
		exit(1);	
	}		
}

sub _getSequenceCoverage{
	my($self)=shift;
	
	if(scalar(@_)==2){
		my $hit=shift;
		my $queryLength=shift;
		return (($hit->hsp_align_len-($hit->hsp_align_len - $hit->hsp_identity)) / $queryLength * 100);
	}
	else{
		print STDERR "incorrect number of arguments sent to getSequenceCoverage!\n";
		exit(1);
	}
}



sub _getCoreAccessoryType{
	my($self)=shift;
	
	if(@_){
		my $result=shift;
		my $returnType;
		my $numberOverSequenceCutoff=0;
		
		#if there is no blast result hit, deal with it
		if(!defined $result->hitHash){
			$returnType = 'accessory';
		}
		else{
			foreach my $hit(keys %{$result->hitHash}){
				my $hitObj = $result->hitHash->{$hit};
				my $sequenceCoverage = $self->_getSequenceCoverage($hitObj,$result->query_len);
				$numberOverSequenceCutoff++ if ($sequenceCoverage >= $self->_percentIdentityCutoff);
			}
			
			if($numberOverSequenceCutoff >= $self->_coreGenomeThreshold){
				$returnType = 'core';
			}
			else{
				$returnType = 'accessory';
			}
		}			
		return $returnType;
	}
	else{
		print STDERR "No item sent to getCoreAccessoryType!\n";
		exit(1);
	}
}



sub _processAccessoryResult{
	my($self)=shift;
	
	my $result =shift;
	my %resultHash;
	my @returnArray;
	
	#add query name
	$returnArray[0]=$result->query_def;
	
	if(defined $result->hitHash){
		foreach my $hit(keys %{$result->hitHash}){
			my $hitObj = $result->hitHash->{$hit};
			my $sequenceCoverage = $self->_getSequenceCoverage($hitObj,$result->query_len);
			
			if($self->_accessoryType eq 'binary'){
				if($sequenceCoverage >= $self->_percentIdentityCutoff){
					$resultHash{$hit}=1;
				}
				else{
					$resultHash{$hit}=0;
				}
			}
			elsif($self->_accessoryType eq 'percent'){
				$resultHash{$hit}=$sequenceCoverage;
			}
			elsif($self->_accessoryType eq 'sequence'){
				$resultHash{$hit}=$hitObj->hsp_hseq;
			}
			else{
				print STDERR "incorrect accessoryType specified!\n";
				exit(1);
			}
		}#end of foreach
	}#end of if
	
	#create the output in correct order
	foreach my $query(keys %{$self->_queryNameOrderHash}){
		my $order = $self->_queryNameOrderHash->{$query};
		
		if(defined $resultHash{$query}){
			$returnArray[$order]=$resultHash{$query}
		}
		else{
			$returnArray[$order]=0;
		}
	}
	
	return (join("\t", @returnArray) . "\n");
}

sub _sendToProcessQueue{
	my($self)=shift;
	
	if(@_){
		my $arrayRef=shift;
		my $resultNumber=shift;
		my @coreOutputBuffer;
		my @accessoryOutputBuffer;
		
		#give each temp process a billion number range to work with
		$resultNumber *= 1000000000;
		
		foreach my $item(@{$arrayRef}){
			$resultNumber++;
			my $type = $self->_getCoreAccessoryType($item);
			
			if($type eq 'core'){
				#get a muscle alignment
				#get tabbed SNPs
				$self->logger->debug("DEBUG:\tProcessing core item $resultNumber");
				push @coreOutputBuffer, $self->_processCoreResult($item,$resultNumber);				
			}
			elsif($type eq 'accessory'){
				#get a tab delimited output in queryName order of either binary, %ID or sequence
				$self->logger->debug("DEBUG:\tProcessing accessory item $resultNumber");
				push @accessoryOutputBuffer, $self->_processAccessoryResult($item);
			}
		}		
		
		if (defined $coreOutputBuffer[0]){
			my $coreTempName = $self->_baseDirectory . 'coreTempFile_' . $resultNumber;	
			my $coreOutFile = IO::File->new('>' . $coreTempName) or die "$!";
			
			foreach my $coreArrayRef(@coreOutputBuffer){
				print $coreOutFile @{$coreArrayRef}; 
			}			
			
			$coreOutFile->close();
			@coreOutputBuffer=();
		}
		
		if (defined $accessoryOutputBuffer[0]){
			my $accessoryTempName = $self->_baseDirectory . 'accessoryTempFile_' . $resultNumber;
			my $accessoryOutFile = IO::File->new('>' . $accessoryTempName) or die "$!";
			print $accessoryOutFile @accessoryOutputBuffer;
			$accessoryOutFile->close();
			@accessoryOutputBuffer=();
		}		
	}
}

sub _processCoreResult{
	my($self)=shift;
	
	my $result=shift;
	my $countNumber=shift;
	
	my @fastaArray;
	my @returnArray;
	
	my %startBpHash;
	foreach my $hit(keys %{$result->hitHash}){
		my $hitObj = $result->hitHash->{$hit};
		$hitObj->setSequence(); #this ensures start <  end bp
		push @fastaArray, ('>' . $hitObj->hit_def . "\n" . $hitObj->hsp_hseq . "\n");
		$startBpHash{$hitObj->hit_def}=$hitObj->hsp_hit_from;	
	}
	
	#create temp files for muscle
	my $tempInFH = File::Temp->new();
	my $tempOutFH = File::Temp->new();
	
	print $tempInFH @fastaArray;
	
	my $muscleCommand = MuscleCmd->new();
	$muscleCommand->setIn($tempInFH->filename);
	$muscleCommand->setOut($tempOutFH->filename);
	$muscleCommand->run();
		
	my @alignedFastaSeqs = $tempOutFH->getlines();
	
	#add SNP information to the return
	my $snpDetective= SNPFinder->new($self->_snpType,$self->_queryNameOrderHash);
	
	push @returnArray, ($snpDetective->findSNPs(\@alignedFastaSeqs,$countNumber,\%startBpHash)); 
	return \@returnArray;
}


sub _combineResultTempFiles{
	my($self)=shift;
	
	if(@_){
		my $numberOfCores=shift;
		my $coreTempFileName = $self->_baseDirectory . 'core_snps_table.txt';
		my $accessoryTempFileName = $self->_baseDirectory . 'accessory_regions_table.txt';
		
		$self->_accessoryTempFile($accessoryTempFileName);
		$self->_coreTempFile($coreTempFileName);
		
		my $forker = Parallel::ForkManager->new($numberOfCores);
		
		my $count=0;
		for(1..2){
			$count++;
			$forker->start and next;
			if($count==1){
				$self->_combineTempFiles('coreTempFile_',$coreTempFileName);
			}
			elsif($count==2){
				$self->_combineTempFiles('accessoryTempFile_',$accessoryTempFileName);
			}			
			$forker->finish;
		}
		$forker->wait_all_children;					
	}
	else{
		print STDERR "FALSE! The sub is a liar\n";
		exit(1);
	}
}

1;
