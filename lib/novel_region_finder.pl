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
use Blast::Annotator;
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
my $mummer = Mummer::MummerGPU->new(
		'baseDirectory'   => $nrf->_baseDirectory
);

#create sub-files for mummer based on the settings in the config file
$mummer->mummersLittleHelper(
	'multiFastaFile'=>$nrf->combinedReferenceFile,
	'bpPerFile'=>$nrf->mummerBpPerSplitFile // undef,
	'numberOfFiles'=>$nrf->mummerNumberOfSplitFiles // undef
);

#allow user specification of number of mummer instances, or default to the number of cores
#allows a fine grain control of memory usage
$mummer->run(
		'queryFile'       => $nrf->combinedQueryFile,
		'mummerDirectory' => $nrf->mummerDirectory,
		'numberOfCores'   => $nrf->mummerNumberOfInstances // $nrf->_numberOfCores
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
if(!defined $nrf->_skipGatherFiles){
	if($nrf->createGraphic eq 'yes'){

			#combine novel regions into single file for proper functioning of SvgDrawer
			my $manipulator = FileInteraction::FileManipulation->new();
			my $singleFile = $nrf->_baseDirectory . 'novelRegions_single.fasta';
			$manipulator->outputFH(IO::File->new('>' . $singleFile)) or die "$!";
			$manipulator->multiFastaToSingleFasta($nrf->outputFileName);

			#run mummer
			#include the base dir in _p
			$logger->info("MUMmer for graphics initiated");
			my $mum = Mummer::MummerGPU->new(
				'baseDirectory'   => $nrf->_baseDirectory
			);

			$mum->run(
					'queryFile'       => $nrf->combinedQueryFile,
					'referenceFile'   => $singleFile,
					'mummerDirectory' => $nrf->mummerDirectory,					
					'numberOfCores'   => $nrf->_numberOfCores,
					'p'=>$nrf->_baseDirectory . 'svg'
			);
			$logger->info("MUMmer for graphics finished");

			#create coords file for svg
			$mum->showCoords(
				'deltaFile'=>$mum->deltaFile,
				'coordsFile'=>$nrf->_baseDirectory . 'graphics.coords'
			);

			#create the visualization
			my $systemLine = "perl $SCRIPT_LOCATION/Visual/svgDrawer.pl " . $mum->coordsFile . ' ' . ' ' . $nrf->outputFileName . ' ' . '1'
				. ' > ' . $nrf->_baseDirectory . 'novelRegionsGraphic.svg';
			$logger->info("Creating graphic with $systemLine");
			system($systemLine);
			$logger->info("Graphic created");
	}

	if($nrf->toAnnotate eq 'yes'){
		$logger->info("Annotation beginning");
		#provide a table of novel region annotations, as well as a newly formatted and annotated novel regions file
		my $annotator = Blast::Annotator->new(
			'inputFile'=>$nrf->outputFileName ,
			'outputFile'=>$nrf->_baseDirectory . '/novel_regions_annotation.txt',
			'blastDirectory'=>$nrf->_blastDirectory,
			'blastDatabase'=>'/home/phac/NCBI/nr',
			'numberOfCores'=>$nrf->_numberOfCores,
			'annotatedFastaFile'=>$nrf->_baseDirectory . 'novelRegions_annotated.fasta'
		);
		$annotator->annotate();
		$logger->info("Annotation complete");
	}
	else{
		$logger->info("No Annotation performed");
	}#end of else
}#end of check for _skipGatherFiles

$nrf->cleanUp();
$nrf->logger->info("Novel regions extracted");



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
