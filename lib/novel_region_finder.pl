#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use File::Basename;
use Modules::Setup::Settings;
use Modules::Setup::PanseqFiles;
use Modules::NovelRegion::NovelIterator;

use Tie::Log4perl;
use Log::Log4perl;

#we need to log in the Settings module
#however, this means we do not know the output directory
#therefore init a basic logger until we have the settings
Log::Log4perl->easy_init();

#get and parse the configuration file for Panseq
my $settings = Modules::Setup::Settings->new($ARGV[0]);

#some programs default messages to STDERR
#we want them logged
#closing STDERR and associating it with Log4perl is done below
close STDERR;
tie *STDERR, "Tie::Log4perl";
Log::Log4perl->init("$FindBin::Bin/log4p.conf");

#get the query/reference files to use in the comparison
my $files = Modules::Setup::PanseqFiles->new(
	'queryDirectory'=>$settings->queryDirectory,
	'referenceDirectory'=>$settings->referenceDirectory
);

my $novelIterator = Modules::NovelRegion::NovelIterator->new(
	'queryFile'=>$files->singleQueryFile($settings->baseDirectory . 'singleQueryFile.fasta'),
	'referenceFile'=>$files->singleReferenceFile($settings->baseDirectory . 'singleReferenceFile.fasta'),
	'panGenomeFile'=>$settings->baseDirectory . 'panGenome.fasta',
	'novelRegionsFile'=>$settings->baseDirectory . 'novelRegions.fasta',
	'settings'=>$settings
);
$novelIterator->run();

# #create visualization if run as novel_region_finder
# #check for _skipGatherFiles to indicate a call from core_accessory.pl
# if(!defined $nrf->_skipGatherFiles){
# 	if($nrf->createGraphic eq 'yes'){

# 			#combine novel regions into single file for proper functioning of SvgDrawer
# 			my $manipulator = FileInteraction::FileManipulation->new();
# 			my $singleFile = $nrf->_baseDirectory . 'novelRegions_single.fasta';
# 			$manipulator->outputFH(IO::File->new('>' . $singleFile)) or die "$!";
# 			$manipulator->multiFastaToSingleFasta($nrf->outputFileName);

# 			#run mummer
# 			#include the base dir in _p
# 			$logger->info("MUMmer for graphics initiated");
# 			my $mum = Mummer::MummerGPU->new(
# 				'baseDirectory'   => $nrf->_baseDirectory
# 			);

# 			$mum->run(
# 					'queryFile'       => $nrf->combinedQueryFile,
# 					'referenceFile'   => $singleFile,
# 					'mummerDirectory' => $nrf->mummerDirectory,					
# 					'numberOfCores'   => $nrf->_numberOfCores,
# 					'p'=>$nrf->_baseDirectory . 'svg'
# 			);
# 			$logger->info("MUMmer for graphics finished");

# 			#create coords file for svg
# 			$mum->showCoords(
# 				'deltaFile'=>$mum->deltaFile,
# 				'coordsFile'=>$nrf->_baseDirectory . 'graphics.coords'
# 			);

# 			#create the visualization
# 			my $systemLine = "perl $SCRIPT_LOCATION/Visual/svgDrawer.pl " . $mum->coordsFile . ' ' . ' ' . $nrf->outputFileName . ' ' . '1'
# 				. ' > ' . $nrf->_baseDirectory . 'novelRegionsGraphic.svg';
# 			$logger->info("Creating graphic with $systemLine");
# 			system($systemLine);
# 			$logger->info("Graphic created");
# 	}

# 	if($nrf->toAnnotate eq 'yes'){
# 		$logger->info("Annotation beginning");
# 		#provide a table of novel region annotations, as well as a newly formatted and annotated novel regions file
# 		my $annotator = Blast::Annotator->new(
# 			'inputFile'=>$nrf->outputFileName ,
# 			'outputFile'=>$nrf->_baseDirectory . '/novel_regions_annotation.txt',
# 			'blastDirectory'=>$nrf->_blastDirectory,
# 			'blastDatabase'=>$nrf->blastDatabase,
# 			'numberOfCores'=>$nrf->_numberOfCores,
# 			'annotatedFastaFile'=>$nrf->_baseDirectory . 'novelRegions_annotated.fasta'
# 		);
# 		$annotator->annotate();
# 		$logger->info("Annotation complete");
# 	}
# 	else{
# 		$logger->info("No Annotation performed");
# 	}#end of else
# 	$nrf->cleanUp();
# }#end of check for _skipGatherFiles

# $nrf->logger->info("Novel region run finished");

# #end



# #HOOKS for log4p.conf
sub nucmer_run_file{
	return ($settings->baseDirectory . 'logs/NucmerRun.pm.log'); 
}

sub panseq_files_file{
	return ($settings->baseDirectory . 'logs/PanseqFiles.pm.log'); 
}

sub novel_region_finder_file{
	return ($settings->baseDirectory . 'logs/NovelRegionFinder.pm.log'); 
}

sub sequence_name_file{
	return ($settings->baseDirectory . 'logs/SequenceName.pm.log'); 
}

sub multi_fasta_sequence_name_file{
	return ($settings->baseDirectory . 'logs/MultiSequenceFastaName.pm.log'); 
}

sub novel_iterator_file{
	return ($settings->baseDirectory . 'logs/NovelIterator.pm.log'); 
}

