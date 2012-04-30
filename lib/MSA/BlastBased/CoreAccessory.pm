#!/usr/bin/perl

package CoreAccessory;

#includes
use strict;
use warnings;
use diagnostics;
use FindBin::libs;
use IO::File;
use FileInteraction::Fasta::SequenceName;
use FileInteraction::SegmentMaker;
use Parallel::ForkManager;
use File::Path qw{make_path remove_tree};
use Blast::FormatBlastDB;
use Blast::BlastIO;
use Blast::BlastParallel;
use Carp;
use FileInteraction::FileManipulation;
use NovelRegion::NovelRegionFinder;
use FileInteraction::Fasta::SequenceRetriever;
use Pipeline::PanseqShared;
use FileInteraction::FlexiblePrinter;
use MSA::BlastBased::CoreAccessoryProcessor;
use TreeBuilding::PhylogenyFileCreator;
use Data::Dumper;
use Logging::Logger;

our @ISA = qw/PanseqShared FlexiblePrinter/;


sub _percentIdentityCutoff{
	my $self=shift;
	$self->{'_CoreAccessory_percentIdentityCutoff'}=shift // return $self->{'_CoreAccessory_percentIdentityCutoff'};
}

sub _coreGenomeThreshold{
	my $self=shift;
	$self->{'_CoreAccessory_coreGenomeThreshold'}=shift // return $self->{'_CoreAccessory_coreGenomeThreshold'};
}

sub _accessoryType{
	my $self=shift;
	$self->{'_CoreAccessory_accessoryType'}=shift // return $self->{'_CoreAccessory_accessoryType'};
}

sub _queryNameOrderHash{
	my $self=shift;
	$self->{'_CoreAccessory_queryNameOrderHash'}=shift // return $self->{'_CoreAccessory_queryNameOrderHash'};
}

sub _snpType{
	my $self=shift;
	$self->{'_CoreAccessory_snpType'}=shift // return $self->{'_CoreAccessory_snpType'};
}

sub _coreInputType{
	my $self=shift;
	$self->{'_CoreAccessory_coreInputType'}=shift // return $self->{'_CoreAccessory_coreInputType'};
}

sub _coreComparisonType{
	my $self=shift;
	$self->{'_CoreAccessory_coreComparisonType'}=shift // return $self->{'_CoreAccessory_coreComparisonType'};
}

sub _segmentCoreInput{
	my $self=shift;
	$self->{'_CoreAccessory_segmentCoreInput'}=shift // return $self->{'_CoreAccessory_segmentCoreInput'};
}

sub _blastType{
	my $self=shift;
	$self->{'_CoreAccessory_blastType'}=shift // return $self->{'_CoreAccessory_blastType'};
}

sub _blastDirectory{
	my $self=shift;
	$self->{'_CoreAccessory_blastDirectory'}=shift // return $self->{'_CoreAccessory_blastDirectory'};
}

sub _muscleExecutable{
	my $self=shift;
	$self->{'_CoreAccessory_muscleDirectory'}=shift // return $self->{'_CoreAccessory_muscleDirectory'};
}

sub _novelSeedName{
	my $self=shift;
	$self->{'_CoreAccessory_novelSeedName'}=shift // return $self->{'_CoreAccessory_novelSeedName'};
}

sub _seedFileName{
	my $self=shift;
	$self->{'_CoreAccessory_seedFileName'}=shift // return $self->{'_CoreAccessory_seedFileName'};
}

sub _notSeedFileName{
	my $self=shift;
	$self->{'_CoreAccessory_notSeedFileName'}=shift // return $self->{'_CoreAccessory_notSeedFileName'};
}


sub _fragmentationSize{
	my $self=shift;
	$self->{'_PanseqShared_fragmentationSize'}=shift // return $self->{'_PanseqShared_fragmentationSize'};
}


sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_coreAccessoryInitialize(@_);
    return $self;
}

#methods
sub _coreAccessoryInitialize{
	my($self)=shift;
	#inheritance
	$self->_panseqSharedInitialize(@_);
	$self->_flexiblePrinterInitialize(@_);
}

#methods

sub _createNovelConfigFile{
	my($self)=shift;
	
	my $fileName = $self->_baseDirectory . 'temp_novel_config.txt';
	my $outFH = IO::File->new('>' . $fileName) or die "$!";
	print $outFH (
		'queryDirectory' . "\t" . $self->queryDirectory . "\n",
		'baseDirectory' . "\t" . $self->_baseDirectory . "\n",
		'numberOfCores' . "\t" . $self->_numberOfCores . "\n",
		'mummerDirectory' . "\t" . $self->mummerDirectory . "\n",
		'minimumNovelRegionSize' . "\t" . $self->minimumNovelRegionSize . "\n",
		'novelRegionFinderMode' . "\t" . 'no_duplicates' . "\n",
		'createGraphic' . "\t" . 'no' . "\n",
		'skipGatherFiles' . "\t" . '1' . "\n",
		'combinedReferenceFile' . "\t" . $self->combinedQueryFile . "\n",
		'combinedQueryFile' . "\t" . $self->_notSeedFileName . "\n"	
	);
		
	$outFH->close();
	return $fileName;
}

sub _createSeedAndNotSeedFiles{
	my($self)=shift;
	
	my $seed = $self->_getSeedFromQueryNameHash();
	
	#create database from queryFiles
	my $retriever = SequenceRetriever->new($self->combinedQueryFile);
	$self->_seedFileName($self->_baseDirectory . 'seedFile_' . $seed->name . '.fasta');
	$self->_notSeedFileName($self->_baseDirectory . 'notSeedFile' . '.fasta');
	
	my $seedFH = IO::File->new('>' . $self->_seedFileName) or die "Cannot open " . $self->_seedFileName . "$!";
	my $notSeedFH = IO::File->new('>' . $self->_notSeedFileName)  or die 'Cannot open ' . $self->_notSeedFileName . "$!";
	
	#create reference file from seed
	foreach my $queryName(keys %{$self->queryNameObjectHash}){
		foreach my $fastaHeader(@{$self->queryNameObjectHash->{$queryName}->arrayOfHeaders}){
			if($self->queryNameObjectHash->{$queryName}->name eq $seed->name){
				print $seedFH '>' . $fastaHeader . "\n" . $retriever->extractRegion($fastaHeader) . "\n";
			}
			else{
				print $notSeedFH '>' . $fastaHeader . "\n" . $retriever->extractRegion($fastaHeader) . "\n";
			}
		}		
	}	
	$seedFH->close();
	$notSeedFH->close();
}

sub runCoreAccessory{
	my($self)=shift;
	
	if(@_){
		my $configFile=shift;
		
		#initialization		
		$self->_validateCoreSettings($self->getSettingsFromConfigurationFile($configFile));
		$self->createDirectories();
		
		my $inputLociFile;
		if($self->_coreInputType eq 'panGenome'){
			$self->getQueryNamesAndCombineAllInputFiles();
			
			$self->logger->debug("DEBUG:\tAfter combining all input files there are: " . (scalar keys %{$self->queryNameObjectHash}));
			$self->logger->info("INFO:\t" . 'Determining non-redundant pan-genome');
			$self->_createSeedAndNotSeedFiles();				
			
			#get no_duplicates novel regions
			my $novelRegionFinder = NovelRegionFinder->new();
									
			#runNovelRegionFinder return file name of novelRegions file
			my $novelRegionFile = $novelRegionFinder->runNovelRegionFinder($self->_createNovelConfigFile);
			
			#combine the novel regions with the seed file for a complete pan-genome
			my $combiner = FileManipulation->new();
			$inputLociFile = $self->_baseDirectory . 'panGenome.fasta';
			my $combinedFH = IO::File->new('>' . $inputLociFile);
			$combiner->outputFilehandle($combinedFH);
			$combiner->vanillaCombineFiles([($self->_seedFileName,$novelRegionFile)]);	
			$combinedFH->close();
		}
		elsif($self->_coreInputType eq 'primers'){
			#set $inputLociFile to fasta file of loci returned from primer search
		}
		elsif($self->_coreInputType =~ /(^\/.+)/){
			$inputLociFile=$1;
			$self->logger->info("INFO:\tUsing " . $inputLociFile . ' as input for Core / Accessory analysis.');
			
			#clean input headers
			my $cleaner = FileManipulation->new();
			
			my $cleanFileName = $self->_baseDirectory . 'cleanedUserProvidedFile.fasta';
			my $cleanFH = IO::File->new('>' . $cleanFileName);
			$cleaner->outputFilehandle($cleanFH);
			$cleaner->combineFilesIntoSingleFile([$inputLociFile],1);
			$cleanFH->close;
			$inputLociFile = $cleanFileName;
			$self->getQueryNamesAndCombineAllInputFiles();
		}
		else{
			print STDERR "incorrect type specified in coreInputType\n";
			exit(1);
		}
		
		#segment the input if required
		if($self->_segmentCoreInput eq 'yes'){
			$self->logger->info("INFO:\tSegmenting the input sequences");
			#create outputFH
			my $segmentedFile = $self->_baseDirectory. 'inputLoci_segments.fasta';
			my $segmentedPanFH = IO::File->new('>' . $segmentedFile) or die "$!";
			
			my $segmenter = SegmentMaker->new($segmentedPanFH);
			$segmenter->segmentTheSequence(
				$inputLociFile,
				$self->_fragmentationSize
			);
			$segmentedPanFH->close();
			$inputLociFile = $segmentedFile;	
		}
		
		#run the comparisons
		if(1){
			my $xmlFiles = $self->_runBlast($inputLociFile,$self->_numberOfCores);
			
			my $cap = CoreAccessoryProcessor->new({
				'baseDirectory'=>$self->_baseDirectory,
				'percentIdentityCutoff'=>$self->_percentIdentityCutoff,
				'queryNameObjectHash'=>$self->queryNameObjectHash,
				'coreGenomeThreshold'=>$self->_coreGenomeThreshold,
				'accessoryType'=>$self->_accessoryType,	
				'snpType'=>$self->_snpType,
				'muscleExecutable'=>$self->_muscleExecutable			
			});
			
			$self->logger->info("INFO:\tBLAST\+ finished. Gathering core \/ accessory information");			
			
			foreach my $file(@{$xmlFiles}){
				$cap->processBlastXML(
					$file,
					$self->_numberOfCores,
					2
				);
			}
			$cap->_combineResultTempFiles($self->_numberOfCores);
		
			#create phylogeny files			
			$self->logger->info("INFO:\tCore \/ accessory information gathered\. Creating phylogeny files\.");
			
			my $forker=Parallel::ForkManager->new($self->_numberOfCores);			
			my $count=0;
			for(1..2){
				$count++;
				$forker->start and next;
				if($count==1){
					$self->_createPhylogenyFiles(
						'phylip',
						$cap->_coreTempFile,
						$self->_baseDirectory . 'core_alignment.phylip',
						$self->_baseDirectory . 'core_alignment_info.txt'	
					);
				}
				elsif($count==2 && $self->_accessoryType eq 'binary'){
					$self->_createPhylogenyFiles(
						'phylip',
						$cap->_accessoryTempFile,
						$self->_baseDirectory . 'accessory_alignment.phylip',
						$self->_baseDirectory . 'accessory_alignment_info.txt'				
					);
				}			
				$forker->finish;
			}
			$forker->wait_all_children;				
		}
		
		$self->logger->info("INFO:\tFinished CoreAccessory Analysis");		
	}
	else{
		print STDERR "incorrect number of arguments sent to runDetermineCoreAccessory\n";
		exit(1);
	}
}

sub _runBlast{ 
	my $self=shift;
	my $inputFile=shift;
	my $numberOfSplits=shift;
	my @blastOutputFileNames;
	
	#create blast db
	#my $blastDB = $self->_createBlastDB();		
	my $blastDB= '/home/phac/all_ecoli/analysis/panseq/no_plasmids_core/output/blastdb';
			
	$self->logger->info("INFO:\tRunning BLAST\+ in seuence");
	
	#blast in parallel works for small datasets, but multi-GB files among multiple processors
	#create an unacceptable memory requirement
	my $splitter = FastaFileSplitter->new();
	$splitter->splitFastaFile($inputFile,$numberOfSplits);
	my $forker = Parallel::ForkManager->new($numberOfSplits);
	my $counter=0;
	
	foreach my $splitFile(@{$splitter->arrayOfSplitFiles}){
		#run blast / set params 
		my $blastFileName = $self->_baseDirectory . $counter . '_blastoutput.xml';
		
		my $blaster = BlastIO->new({
			'blastDirectory'=>$self->_blastDirectory,
			'type'=>$self->_blastType,
			'db'=>$blastDB,
			'outfmt'=>'5',
			'evalue'=>'0.00001',
			'word_size'=>20,
			'num_threads'=>$self->_numberOfCores
		});
		
		
		push @blastOutputFileNames, $blastFileName;
		$counter++;
		
		$forker->start and next;
			$blaster->runBlastn(
				'query'=>$splitFile,
				'out'=>$blastFileName
			);
		$forker->finish();	
		
	}
	$forker->wait_all_children();
	return \@blastOutputFileNames;
}

sub _createPhylogenyFiles{
	my($self)=shift;
	
	if(scalar(@_)==4){
		my $type=shift;
		my $tabFile=shift;
		my $phylogenyOutputName = shift;
		my $phylogenyOutputInfoName = shift;
		
		$self->logger->info("INFO:\tCreating phylogeny files $phylogenyOutputName and $phylogenyOutputInfoName");
		
		#check to see if file exists and is non-zero in size
		if(-s $tabFile){
		
			#create core phylogeny file based on snps
			my $phyllo= PhylogenyFileCreator->new();
			my $phylogenyFH = IO::File->new('>' . $phylogenyOutputName) or die "$!";
			my $phylogenyInfoFH = IO::File->new('>' . $phylogenyOutputInfoName) or die "$!";
			
			$phyllo->outputFilehandle($phylogenyFH);
			$phyllo->phylipInfoFH($phylogenyInfoFH); #defaults to STDOUT
			
			#requires <type>,<input file>, <header flag 0/1> (optional. default=0)
			$phyllo->createFile(
				$type,
				$tabFile,
				1
			);
			$phylogenyFH->close();
			$phylogenyInfoFH->close();
		}
	}
		else{
		print STDERR "Wrong number of arguments sent to createPhylogenyFiles!\n";
		exit(1);
	}
}



sub _createBlastDB{
	my($self)=shift;	
	
	my $dbMaker = FormatBlastDB->new($self->_blastDirectory);
	
	my $blastDB = $self->_baseDirectory . 'blastdb';	
	
	$self->logger->info("INFO:\tCreating the blast database as $blastDB");
	
	$dbMaker->runMakeBlastDb(
		'dbtype'=>'nucl',
		'in'=>$self->combinedQueryFile,		
		'out'=>$blastDB,
		'title'=>'blastdb'	
	);
	return $blastDB;
}



sub _validateCoreSettings {
	my ($self) = shift;

	my $validator = $self->_validator;

	if (@_) {
		my $settingsHashRef = shift;

		foreach my $setting ( keys %{$settingsHashRef} ) {
			my $value = $settingsHashRef->{$setting};

			$self->_segmentCoreInput( $validator->yesOrNo($value) )       if $setting eq 'segmentCoreInput';
			$self->_blastDirectory( $validator->isADirectory($value) )    if $setting eq 'blastDirectory';
			$self->_muscleExecutable( $validator->doesFileExist($value) ) if $setting eq 'muscleExecutable';

			#unique checks
			$self->_percentIdentityCutoff( $validator->isAValidPercentID($value)*100 ) if $setting eq 'percentIdentityCutoff';
			$self->_coreGenomeThreshold( $validator->isAnIntGreaterThan( $value, -1 ) ) if $setting eq 'coreGenomeThreshold';
			$self->_novelRegionFinderMode( $validator->novelRegionFinderModeCheck($value) ) if $setting eq 'novelRegionFinderMode';
			$self->_accessoryType( $validator->accessoryTypeCheck($value) )                 if $setting eq 'accessoryType';
			$self->_coreInputType($value)                                                   if $setting eq 'coreInputType';
			#set snpType
			$self->_coreComparisonType( $validator->coreComparisonTypeCheck($value) )       if $setting eq 'coreComparisonType';
			if ( $setting eq 'blastType' ) {
				$self->_blastType( $validator->blastTypeCheck($value) );
				if ( $value eq 'blastn' ) {
					$self->_snpType('nucleotide');
				}
				elsif ( $value eq 'tblastn' ) {
					$self->_snpType('protein');
				}
			}
			$self->_fragmentationSize( $validator->isAnInt($value) ) if $setting eq 'fragmentationSize';
		}
	}
}

sub _coreGenomeThresholdCheck{
	my($self)=shift;
	
	if(@_){
		my $thresh=shift;
		
		if(($thresh eq 'all') || ($self->isAnInt($thresh))){
			return $thresh;
		}
		else{
			print STDERR "$thresh is not a valid coreGenomeTheshold value!\n",
				"Value needs to be 'all' or a positive integer.\n";
			exit(1);
		}
	}
	else{
		print STDERR "nothing sent to coreGenomeThresholdCheck\n";
		exit(1);
	}
}

sub _accessoryTypeCheck{
	my($self)=shift;
	
	if(@_){
		my $type=shift;
		
		if(($type eq 'binary') || ($type eq 'percent') || ($type eq 'sequence')){
			return $type;
		}
		else{
			print STDERR "$type is not a valid accessoryType!\n",
				"Valid options are binary, percent and sequence!\n";
			exit(1);
		}
	}
	else{
		print STDERR "nothing sent to accessoryTypeCheck!\n";
		exit(1);
	}
}

sub _coreComparisonTypeCheck{
	my($self)=shift;
	
	if(@_){
		my $type=shift;
		
		if(($type eq 'blast') || ($type eq 'nucmer')){
			return $type;
		}
		else{
			print STDERR "$type is not a valid coreComparisonType!\n",
				"Valid options are blast and nucmer!\n";
			exit(1);
		}
	}
	else{
		print STDERR "nothing sent to accessoryTypeCheck!\n";
		exit(1);
	}
}


sub _blastTypeCheck{
	my($self)=shift;
	
	if(@_){
		my $type=shift;
		
		if($type eq 'blastn'){
			$self->_snpType('nucleotide');
			return $type;
		}
		elsif($type eq 'tblastn'){
			$self->_snpType('protein');
			return $type;
		}
		else{
			$self->logger->fatal("$type is an invalid entry for blastType in the configuration file!\n
				Currently supported values are blastn and tblastn");
			exit(1);
		}
	}
	else{
		print STDERR "Nothing sent to blastTypeCheck!\n";
		exit(1);
	}
}

sub _getSeedFromQueryNameHash{
	my($self)=shift;
		
	my $seed;
	foreach my $queryName(keys %{$self->queryNameObjectHash}){
		$seed=$queryName;
		my $numberOfFastaSequences = scalar(@{$self->queryNameObjectHash->{$queryName}->arrayOfHeaders});
		
		#if single fasta header, means a closed genome
		if($numberOfFastaSequences ==1){
			last;
		}
	}	
	$self->logger->info("INFO:\t" . $self->queryNameObjectHash->{$seed}->name . ' selected as seed sequence.');
		
	return $self->queryNameObjectHash->{$seed};
}

1;
