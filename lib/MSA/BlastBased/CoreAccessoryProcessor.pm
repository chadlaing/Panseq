#!/usr/bin/perl

package MSA::BlastBased::CoreAccessoryProcessor;

use strict;
use warnings;
use diagnostics;
use Carp;
use FindBin;
use lib "$FindBin::Bin";
use MSA::BlastBased::BlastResultFactory;
use MSA::BlastBased::SNPFinder;
use Muscle::MuscleCmd;
use FileInteraction::FileManipulation;
use File::Temp;

use parent 'MSA::BlastBased::CoreAccessory';

sub new{
	my $class = shift;
	my $self= $class->SUPER::new(@_); #this calls CoreAccessoryProcessor::_initialize, so need to call the SUPER::_initialize
	return $self;
}

#class variables
sub _queryNameOrderHash {
	my $self = shift;
	$self->{'_CoreAccessoryProcessor_queryNameOrderHash'} = shift // return $self->{'_CoreAccessoryProcessor_queryNameOrderHash'};
}

sub _accessoryTempFile {
	my $self = shift;
	$self->{'_CoreAccessoryProcessor_accessoryTempFile'} = shift // return $self->{'_CoreAccessoryProcessor_accessoryTempFile'};
}

sub _coreTempFile {
	my $self = shift;
	$self->{'_CoreAccessoryProcessor_coreTempFile'} = shift // return $self->{'_CoreAccessoryProcessor_coreTempFile'};
}


#methods
sub _initialize {
	my $self = shift;

	#inheritance
	$self->SUPER::_initialize(@_);

	#logging
	$self->logger(Log::Log4perl->get_logger());

	#set values
	my $paramsRef = shift;
	$self->_baseDirectory( $paramsRef->{'baseDirectory'} )                 // confess('baseDirectory required in coreAccessoryProcessor');
	$self->_percentIdentityCutoff( $paramsRef->{'percentIdentityCutoff'} ) // confess('percentIdentityCutoff required in coreAccessoryProcessor');
	$self->queryNameObjectHash( $paramsRef->{'queryNameObjectHash'} )      // confess('queryNameObjectHash required in coreAccessoryProcessor');
	$self->_coreGenomeThreshold( $paramsRef->{'coreGenomeThreshold'} )     // confess('coreGenomeThreshold required in coreAccessoryProcessor');
	$self->_accessoryType( $paramsRef->{'accessoryType'} )                 // confess('accessoryType required in coreAccessoryProcessor');
	$self->_snpType( $paramsRef->{'snpType'} )                             // confess('snpType required in coreAccessoryProcessor');
	$self->_muscleExecutable( $paramsRef->{'muscleExecutable'} )           // confess('muscleExecutable required in coreAccessoryProcessor');

	#set up queryNameOrderHash
	$self->_getQueryNameOrder();


}

sub _getQueryNameOrder {
	my ($self) = shift;

	my %returnHash;
	my $order = 1;

	$self->logger->debug( "Ref type of queryNameObjectHash" . ref( $self->queryNameObjectHash ) );
	$self->logger->debug( "Scalar of queryNameObjectHash" . scalar keys %{ $self->queryNameObjectHash } );

	foreach my $name ( keys %{ $self->queryNameObjectHash } ) {
		$returnHash{$name} = $order;
		$order++;
	}

	$self->_queryNameOrderHash( \%returnHash );
}

sub _combineTempFiles {
	my ($self) = shift;

	if ( scalar(@_) == 2 ) {
		my $textToMatch    = shift;
		my $outputFileName = shift;

		my @filesToCombine;
		my $fm       = FileInteraction::FileManipulation->new();
		my $namesRef = $fm->getFileNamesFromDirectory( $self->_baseDirectory );

		foreach my $name ( @{$namesRef} ) {
			if ( $name =~ /$textToMatch\d+/ ) {
				push @filesToCombine, ($name);
			}
		}

		if ( scalar(@filesToCombine) > 0 ) {
			my $oFile = IO::File->new( '>' . $outputFileName ) or die "$!";

			#print out header names on combined files (tab-delimited names)
			my @outputArray;
			foreach my $name ( keys %{ $self->_queryNameOrderHash } ) {
				$outputArray[ $self->_queryNameOrderHash->{$name} ] = $name;
			}

			for ( my $i = 1 ; $i < scalar(@outputArray) ; $i++ ) {
				print $oFile "\t" . $outputArray[$i];
			}
			print $oFile "\n";

			$fm->outputFH($oFile);
			$fm->vanillaCombineFiles( \@filesToCombine,1 ); #1 for destroy
			$oFile->close();
		}
	}
	else {
		$self->logger->fatal("Wrong number of arguments to combineTempFiles!");
		exit(1);
	}
}

sub processBlastXML {
	my ($self) = shift;

	if ( scalar(@_) == 3 ) {
		my $blastXMLFile    = shift;
		my $numberOfCores   = shift;
		my $resultArraySize = shift;

		$self->logger->info("Processing Blast XML file $blastXMLFile with $numberOfCores cores and resultArraySize of $resultArraySize");

		my $forker        = Parallel::ForkManager->new($numberOfCores);
		my $currentResult = 0;
		my @resultArray;

		#create filehandle
		my $xmlFH = IO::File->new( '<' . $blastXMLFile ) or die "$!";
		my $xmler = MSA::BlastBased::BlastResultFactory->new($xmlFH);

		while ( my $result = $xmler->nextResult ) {
			push @resultArray, $result;

			if ( scalar(@resultArray) == $resultArraySize ) {
				$self->logger->debug("Starting fork. Result $currentResult");
				$forker->start and next;
				$self->_sendToProcessQueue( \@resultArray, $currentResult, $resultArraySize,$blastXMLFile );
				$forker->finish;
			}
		}
		continue {
			$currentResult++;
			if ( scalar(@resultArray) == $resultArraySize ) {
				$self->logger->debug("Result array emptied");
				@resultArray = ();
			}
		}
		#$account for stored but unprocessed items
		if(scalar(@resultArray) > 0){
			$self->_sendToProcessQueue(\@resultArray, $currentResult, $resultArraySize,$blastXMLFile );
		}
		
		$forker->wait_all_children();
		$xmlFH->close();

	}    #end if
	else {
		$self->logger->fatal("please specify a BLAST xml output file for processing!");
		exit(1);
	}
}

sub _getSequenceCoverage {
	my ($self) = shift;

	if ( scalar(@_) == 2 ) {
		my $hit         = shift;
		my $queryLength = shift;
		return ( ( $hit->hsp_align_len - ( $hit->hsp_align_len - $hit->hsp_identity ) ) / $queryLength * 100 );
	}
	else {
		$self->logger->fatal("incorrect number of arguments sent to getSequenceCoverage!");
		exit(1);
	}
}

sub _getCoreAccessoryType {
	my ($self) = shift;

	if (@_) {
		my $result = shift;
		my $returnType;
		my $numberOverSequenceCutoff = 0;

		#if there is no blast result hit, deal with it
		if ( !defined $result->hitHash ) {
			$self->logger->debug("No blast result hit for result");
			$returnType = 'accessory';
		}
		else {
			foreach my $hit ( keys %{ $result->hitHash } ) {
				my $hitObj = $result->hitHash->{$hit};
				my $sequenceCoverage = $self->_getSequenceCoverage( $hitObj, $result->query_len );
				$numberOverSequenceCutoff++ if ( $sequenceCoverage >= $self->_percentIdentityCutoff );
			}

			if ( $numberOverSequenceCutoff >= $self->_coreGenomeThreshold ) {
				$returnType = 'core';
			}
			else {
				$returnType = 'accessory';
			}
		}
		$self->logger->debug("TYPE: $returnType number over cutoff: $numberOverSequenceCutoff");
		return $returnType;
	}
	else {
		$self->logger->fatal("No item sent to getCoreAccessoryType!");
		exit(1);
	}
}

sub _processAccessoryResult {
	my ($self) = shift;

	my $result = shift;
	my %resultHash;
	my @returnArray;

	#add query name
	$returnArray[0] = $result->query_def;

	if ( defined $result->hitHash ) {
		foreach my $hit ( keys %{ $result->hitHash } ) {
			my $hitObj = $result->hitHash->{$hit};
			my $sequenceCoverage = $self->_getSequenceCoverage( $hitObj, $result->query_len );

			if ( $self->_accessoryType eq 'binary' ) {
				if ( $sequenceCoverage >= $self->_percentIdentityCutoff ) {
					$resultHash{$hit} = 1;
				}
				else {
					$resultHash{$hit} = 0;
				}
			}
			elsif ( $self->_accessoryType eq 'percent' ) {
				$resultHash{$hit} = $sequenceCoverage;
			}
			elsif ( $self->_accessoryType eq 'sequence' ) {
				$resultHash{$hit} = $hitObj->hsp_hseq;
			}
			else {
				$self->logger->fatal("incorrect accessoryType specified!");
				exit(1);
			}
		}    #end of foreach
	}    #end of if

	#create the output in correct order
	foreach my $query ( keys %{ $self->_queryNameOrderHash } ) {
		my $order = $self->_queryNameOrderHash->{$query};

		if ( defined $resultHash{$query} ) {
			$returnArray[$order] = $resultHash{$query};
		}
		else {
			$returnArray[$order] = 0;
		}
	}

	return ( join( "\t", @returnArray ) . "\n" );
}

sub _sendToProcessQueue {
	my ($self) = shift;

	if (@_) {
		my $arrayRef        = shift;
		my $resultNumber    = shift;
		my $resultArraySize = shift;
		my $xmlFileName 	=shift;
		my @coreOutputBuffer;
		my @accessoryOutputBuffer;

		#give each temp process a billion number range to work with
		$resultNumber *= 1000000000;

		my $purgeNumber = 0;
		foreach my $item ( @{$arrayRef} ) {
			$purgeNumber++;
			$resultNumber++;
			my $type = $self->_getCoreAccessoryType($item);
			
			if ( $type eq 'core' ) {

				#get a muscle alignment
				#get tabbed SNPs
				$self->logger->debug("Processing core item $resultNumber"); 
				my $resultArrayRef = $self->_processCoreResult( $item, $resultNumber );
				if(scalar @{$resultArrayRef} > 0){
					push @coreOutputBuffer, $resultArrayRef;
				}
				
			}
			elsif ( $type eq 'accessory' ) {

				#get a tab delimited output in queryName order of either binary, %ID or sequence
				$self->logger->debug("Processing accessory item $resultNumber");
				my $result = $self->_processAccessoryResult( $item);
				push @accessoryOutputBuffer, $result;
				
			}

			if ( $purgeNumber == $resultArraySize ) {
				$self->_purgeBuffer( \@coreOutputBuffer, \@accessoryOutputBuffer, $resultNumber,$xmlFileName );
				$purgeNumber           = 0;
				@accessoryOutputBuffer = ();
				@coreOutputBuffer      = ();
			}
		}
	}
}

sub _purgeBuffer {
	my $self         = shift;
	
	my $cRef         = shift;
	my $aRef         = shift;
	my $resultNumber = shift;
	my $xmlFileName	=shift;
	
	#remove non-word chars from filename
	$xmlFileName =~ s/\W//g;

	if ( defined $cRef->[0]) {
		my $coreTempName = $self->_baseDirectory . $xmlFileName . '_coreTempFile_' . $resultNumber;
		my $coreOutFile = IO::File->new( '>' . $coreTempName ) or die "$!";

		foreach my $coreArrayRef ( @{$cRef} ) {
			foreach my $coreLineArrayRef(@{$coreArrayRef}){
				print $coreOutFile @{$coreLineArrayRef};
			}
		}
		$coreOutFile->close();
	}
	else {
		$self->logger->debug("core ref is empty at result $resultNumber");
	}

	if ( defined $aRef->[0] ) {
		my $accessoryTempName = $self->_baseDirectory . $xmlFileName .'_accessoryTempFile_' . $resultNumber;
		my $accessoryOutFile = IO::File->new( '>' . $accessoryTempName ) or die "$!";
		print $accessoryOutFile @{$aRef};
		$accessoryOutFile->close();
	}
	else {
		$self->logger->debug("accessory ref is empty at result $resultNumber");
	}
}

sub _processCoreResult {
	my ($self) = shift;

	my $result      = shift;
	my $countNumber = shift;

	my @fastaArray;
	my @returnArray;

	my %startBpHash;
	foreach my $hit ( keys %{ $result->hitHash } ) {
		my $hitObj = $result->hitHash->{$hit};
		$hitObj->setSequence();    #this ensures start <  end bp
		push @fastaArray, ( '>' . $hitObj->hit_def . "\n" . $hitObj->hsp_hseq . "\n" );
		$startBpHash{ $hitObj->hit_def } = $hitObj->hsp_hit_from;
	}

	#create temp files for muscle
	my $tempInFH  = File::Temp->new();
	my $tempOutFH = File::Temp->new();

	print $tempInFH @fastaArray;

	my $muscleCommand = Muscle::MuscleCmd->new( $self->_muscleExecutable );
	$muscleCommand->setIn( $tempInFH->filename );
	$muscleCommand->setOut( $tempOutFH->filename );
	$muscleCommand->run();

	my @alignedFastaSeqs = $tempOutFH->getlines();

	#add SNP information to the return
	my $snpDetective = MSA::BlastBased::SNPFinder->new( $self->_snpType, $self->_queryNameOrderHash );
	
	my $snpDetectiveResultArrayRef = $snpDetective->findSNPs( \@alignedFastaSeqs, $countNumber, \%startBpHash );
	
	#check that there are actually SNPs, otherwise return undef
	if(scalar(@{$snpDetectiveResultArrayRef})>0){ 
		push @returnArray, $snpDetectiveResultArrayRef;
		$self->logger->debug("Returning " . @{$snpDetectiveResultArrayRef} . ' lines');
	}
	else{
		$self->logger->debug("SNPFinder empty");
	}
	return \@returnArray;
}

sub _combineResultTempFiles {
	my ($self) = shift;

	if (@_) {
		my $numberOfCores         = shift;
		my $coreTempFileName      = $self->_baseDirectory . 'core_snps_table.txt';
		my $accessoryTempFileName = $self->_baseDirectory . 'accessory_regions_table.txt';

		$self->_accessoryTempFile($accessoryTempFileName);
		$self->_coreTempFile($coreTempFileName);

		my $forker = Parallel::ForkManager->new($numberOfCores);

		my $count = 0;
		for ( 1 .. 2 ) {
			$count++;
			$forker->start and next;
			if ( $count == 1 ) {
				$self->_combineTempFiles( 'coreTempFile_', $coreTempFileName );
			}
			elsif ( $count == 2 ) {
				$self->_combineTempFiles( 'accessoryTempFile_', $accessoryTempFileName );
			}
			$forker->finish;
		}
		$forker->wait_all_children;
	}
	else {
		$self->logger->fatal("FALSE! The sub is a liar");
		exit(1);
	}
}

1;
