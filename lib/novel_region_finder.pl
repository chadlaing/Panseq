#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use FindBin;
use lib "$FindBin::Bin";
use File::Path qw{make_path remove_tree};
use File::Basename;
use NovelRegion::NovelRegionFinder;
use FileInteraction::Fasta::SequenceRetriever;
use Mummer::MummerGPU;
use FileInteraction::FileManipulation;
use IO::File;
#use Visual::Visualization;
#use Blast::Annotator;
use Log::Log4perl qw/get_logger/;
use Tie::Log4perl;

#get script location via File::Basename
my $SCRIPT_LOCATION = dirname(__FILE__);

#usage: novelRegionFinder.pl <config file>
#verify config file parameter
my $configFile;
if(defined $ARGV[0]){
	$configFile=$ARGV[0];
}
else{
	print STDERR "Please specify a configuration file!\n";
	exit(1);
}

my $nrf = NovelRegion::NovelRegionFinder->new();
#initialization
$nrf->validateNovelSettings($nrf->getSettingsFromConfigurationFile($configFile));
my $BASE_DIR=$nrf->_baseDirectory;

unless(defined $nrf->_skipGatherFiles) {
	if ( defined $nrf->_baseDirectory ) {		
		#uses File::Path
		remove_tree( $nrf->_baseDirectory );
		make_path( $nrf->_baseDirectory );
		make_path( $nrf->_baseDirectory . '/logs/');
	}
	else {
		print STDERR "_baseDirectory undefined!\n";
		exit(1);
	}
}

#MUMmmer defaults its messages to STDERR
#we want them logged
#closing STDERR and associating it with Log4perl is done below
#the logger.MummerGPU section in log4p.conf logs this output to a file

close STDERR;
tie *STDERR, "Tie::Log4perl";
Log::Log4perl->init("$SCRIPT_LOCATION/Logging/log4p.conf");
my $logger = Log::Log4perl->get_logger();

$logger->info("Starting novel_region_finder with configFile: $configFile");
$logger->info("Starting novel_region_finder with BASE_DIR: $BASE_DIR");

#this allows the NovelRegionFinder to be called without having to recombine files
#and create directories
if ( $nrf->_skipGatherFiles ) {
	my $fileManipulator = FileInteraction::FileManipulation->new();

	$logger->info( "Skip gathering files. Using previously created combined query file: " . $nrf->combinedQueryFile );
	$nrf->queryNameObjectHash( $nrf->getSequenceNamesAsHashRef( $fileManipulator->getFastaHeadersFromFile( $nrf->combinedQueryFile ) ) );
}
else {
	$nrf->getQueryNamesAndCombineAllInputFiles();
	$logger->info("Created combined query file");
}

#set up the mummer run
my $mummer = Mummer::MummerGPU->new();
$mummer->run(
		'queryFile'       => $nrf->combinedQueryFile,
		'referenceFile'   => $nrf->combinedReferenceFile,
		'mummerDirectory' => $nrf->mummerDirectory,
		'baseDirectory'   => $nrf->_baseDirectory,
		'numberOfCores'   => $nrf->_numberOfCores
);
$logger->info("Nucmer finished. Gathering novel regions.");

$nrf->findNovelRegions( { 'deltaFile' => $mummer->deltaFile, } );
$logger->info("Novel regions gathered. Extracting.");


#open output filehandle
$nrf->outputFileName($nrf->_baseDirectory . 'novelRegions.fasta');
#order: <novel hash ref>, minimumNovelRegionSize, <output FH if different from STDOUT>
my $gr = FileInteraction::Fasta::SequenceRetriever->new( 
	'inputFile'=>$nrf->combinedQueryFile,
	'outputFile'=>$nrf->outputFileName 
);
$gr->extractAndPrintRegionsFromHash(
		'hashRef' => $nrf->novelRegionsHashRef,
		'cutoffSize' => $nrf->minimumNovelRegionSize
);
$logger->info("Novel regions extracted.");

#create visualization if run as novel_region_finder
#check for _skipGatherFiles to indicate a call from core_accessory.pl
# if(!defined $nrf->_skipGatherFiles){
# 	if($nrf->createGraphic eq 'yes'){
# 		my $visFH = IO::File->new('>' . $nrf->_baseDirectory . 'novelRegions.svg') or die "$!";
# 		my $vis = Visual::Visualization->new();
# 		$visFH->print($vis->run($mummer->deltaFile));
# 		$visFH->close;	
# 	}

# 	#provide a table of novel region annotations, as well as a newly formatted and annotated novel regions file
# 	my $annotator = Blast::Annotator->new(
# 		'inputFile'=>$nrf->outputFileName ,
# 		'outputFile'=>$nrf->_baseDirectory . '/novel_regions_annotation.txt',
# 		'blastDirectory'=>'/home/phac/ncbi-blast-2.2.26+/bin/',
# 		'blastDatabase'=>'/home/phac/workspace/Panseq_dev/Panseq2/NCBI_DB/NR',
# 		'numberOfCores'=>'20',
# 		'annotatedFastaFile'=>$nrf->_baseDirectory . 'novelRegions_annotated.fasta'
# 	);
# 	$annotator->annotate();
# }





#HOOKS for log4p.conf
sub novel_region_finder_file{
	return ($BASE_DIR . 'logs/NovelRegionFinder.pm.log'); 
}

sub mummer_gpu_file{
	return ($BASE_DIR . 'logs/MummerGPU.pm.log'); 
}

sub file_manipulation_file{
	return ($BASE_DIR . 'logs/FileManipulation.pm.log'); 
}
