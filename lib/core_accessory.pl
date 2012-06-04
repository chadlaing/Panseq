#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use FindBin::libs;
use MSA::BlastBased::CoreAccessory;
use Log::Log4perl;
use Tie::Log4perl;
use File::Basename;
use File::Slurp;
use File::Path qw{make_path remove_tree};

#get script location via File::Basename
my $SCRIPT_LOCATION = dirname(__FILE__);

=head1
Takes in a single directory name and configuration file as input
Outputs all files to the directory specified in the configuration file.
=cut

my $configFile =  $ARGV[0];
unless(-e $configFile){
	print STDERR "Configuration file not specified or does not exist!";
	exit(1);
}

#run
my $ca = CoreAccessory->new();

#initialization		
$ca->_validateCoreSettings($ca->getSettingsFromConfigurationFile($configFile));
my $BASE_DIR=$ca->_baseDirectory;
if ( defined $ca->_baseDirectory ) {		
		#uses File::Path
		remove_tree( $ca->_baseDirectory );
		make_path( $ca->_baseDirectory );
		make_path( $ca->_baseDirectory . '/logs/');
	}
	else {
		print STDERR "_baseDirectory undefined!\n";
		exit(1);
	}

#closing STDERR and associating it with Log4perl is done below
#the logger.MummerGPU section in log4p.conf logs this output to a file
close STDERR;
tie *STDERR, "Tie::Log4perl";
Log::Log4perl->init("$SCRIPT_LOCATION/Logging/core_accessory_log4p.conf");
my $logger = Log::Log4perl->get_logger();
		
my $inputLociFile;
if($ca->_coreInputType eq 'panGenome'){
	$ca->getQueryNamesAndCombineAllInputFiles();
		
	$ca->logger->debug("After combining all input files there are: " . (scalar keys %{$ca->queryNameObjectHash}));
	$ca->logger->info('Determining non-redundant pan-genome');
	$ca->_createSeedAndNotSeedFiles();				
			
	#get no_duplicates novel regions
												
	#runNovelRegionFinder return file name of novelRegions file
	system('perl ' . $SCRIPT_LOCATION. '/novel_region_finder.pl ' . $ca->_createNovelConfigFile);
	my $novelRegionFile = $ca->_baseDirectory . 'novelRegions.fasta';
			
	#combine the novel regions with the seed file for a complete pan-genome
	my $combiner = FileManipulation->new();
	$inputLociFile = $ca->_baseDirectory . 'panGenome.fasta';
	my $combinedFH = IO::File->new('>' . $inputLociFile);
	$combiner->outputFilehandle($combinedFH);
	$combiner->vanillaCombineFiles([($ca->_seedFileName,$novelRegionFile)]);	
	$combinedFH->close();
}
elsif($ca->_coreInputType eq 'primers'){
	#set $inputLociFile to fasta file of loci returned from primer search
}
elsif($ca->_coreInputType =~ /(^\/.+)/){
	$inputLociFile=$1;
	$ca->logger->info("INFO:\tUsing " . $inputLociFile . ' as input for Core / Accessory analysis.');
		
	#clean input headers
	my $cleaner = FileManipulation->new();
			
	my $cleanFileName = $ca->_baseDirectory . 'cleanedUserProvidedFile.fasta';
	my $cleanFH = IO::File->new('>' . $cleanFileName) or die "Cannot open $cleanFileName $!";
	$cleaner->outputFilehandle($cleanFH);
	$cleaner->combineFilesIntoSingleFile([$inputLociFile],1);
	$cleanFH->close;
	$inputLociFile = $cleanFileName;
	$ca->queryNameObjectHash( $ca->getSequenceNamesAsHashRef( $cleaner->getFastaHeadersFromFile( $inputLociFile ) ) );
}
else{
	$ca->logger->fatal("incorrect type specified in coreInputType");
	exit(1);
}
		
#segment the input if required
if($ca->_segmentCoreInput eq 'yes'){
	$ca->logger->info("Segmenting the input sequences");
	#create outputFH
	my $segmentedFile = $ca->_baseDirectory. 'inputLoci_segments.fasta';
	my $segmentedPanFH = IO::File->new('>' . $segmentedFile) or die "$!";
		
	my $segmenter = SegmentMaker->new($segmentedPanFH);
	$segmenter->segmentTheSequence(
		$inputLociFile,
		$ca->_fragmentationSize
	);
	$segmentedPanFH->close();
	$inputLociFile = $segmentedFile;	
}
		
#run the comparisons
my $xmlFiles = $ca->_runBlast($inputLociFile,$ca->_numberOfCores);
my $cap = CoreAccessoryProcessor->new({
	'baseDirectory'=>$ca->_baseDirectory,
	'percentIdentityCutoff'=>$ca->_percentIdentityCutoff,
	'queryNameObjectHash'=>$ca->queryNameObjectHash,
	'coreGenomeThreshold'=>$ca->_coreGenomeThreshold,
	'accessoryType'=>$ca->_accessoryType,	
	'snpType'=>$ca->_snpType,
	'muscleExecutable'=>$ca->_muscleExecutable			
});
			
$ca->logger->info("BLAST\+ finished. Gathering core \/ accessory information");			
			
foreach my $file(@{$xmlFiles}){
	$cap->processBlastXML(
		$file,
		$ca->_numberOfCores,
		2
	);
}
$cap->_combineResultTempFiles($ca->_numberOfCores);
		
#create phylogeny files			
$ca->logger->info("Core \/ accessory information gathered\. Creating phylogeny files\.");
my $forker=Parallel::ForkManager->new($ca->_numberOfCores);			
my $count=0;
for(1..2){
	$count++;
	$forker->start and next;
	if($count==1){
		$ca->_createPhylogenyFiles(
		'phylip',
		$cap->_coreTempFile,
		$ca->_baseDirectory . 'core_alignment.phylip',
		$ca->_baseDirectory . 'core_alignment_info.txt'	
		);
	}
	elsif($count==2 && $ca->_accessoryType eq 'binary'){
		$ca->_createPhylogenyFiles(
		'phylip',
		$cap->_accessoryTempFile,
		$ca->_baseDirectory . 'accessory_alignment.phylip',
		$ca->_baseDirectory . 'accessory_alignment_info.txt'				
		);
	}			
	$forker->finish;
}
$forker->wait_all_children;				
$ca->logger->info("Finished CoreAccessory Analysis");		



#HOOKS for log4p.conf
sub fasta_file_splitter_file{
	return ($BASE_DIR . 'logs/FastaFileSplitter.pm.log'); 
}

sub blast_parallel_file{
	return ($BASE_DIR . 'logs/BlastParallel.pm.log'); 
}

sub segment_maker_file{
	return ($BASE_DIR . 'logs/BlastParallel.pm.log'); 
}

sub format_blast_db_file{
	return ($BASE_DIR . 'logs/FormatBlastDB.pm.log'); 
}

sub file_manipulation_file{
	return ($BASE_DIR . 'logs/FileManipulation.pm.log'); 
}

sub phylogeny_file_creator_file{
	return ($BASE_DIR . 'logs/PhylogenyFileCreator.pm.log'); 
}

sub blast_result_factory_file{
	return ($BASE_DIR . 'logs/BlastResultFactory.pm.log')
}

sub core_accessory_processor_file{
	return ($BASE_DIR . 'logs/CoreAccessoryProcessor.pm.log')
}

sub blast_result_object_file{
	return ($BASE_DIR . 'logs/BlastResultObject.pm.log')
}

sub blast_hit_object_file{
	return ($BASE_DIR . 'logs/BlastHitObject.pm.log')
}

sub snp_finder_file{
	return ($BASE_DIR . 'logs/SNPFinder.pm.log')
}

sub muscle_cmd_file{
	return ($BASE_DIR . 'logs/MuscleCmd.pm.log')
}

