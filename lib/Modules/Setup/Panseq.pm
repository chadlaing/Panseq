#!/usr/bin/env perl

=pod

=head1 NAME

Modules::Setup::Panseq 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: 

=head1 AUTHOR

Your name (yourname@email.com)

=head2 Methods

=cut

package Modules::Setup::Panseq;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use File::Basename;
use Modules::Setup::PanseqFiles;
use Modules::NovelRegion::NovelIterator;
use Modules::Alignment::BlastRun;
use Modules::Alignment::MakeBlastDB;
use Modules::Fasta::SegmentMaker;
use Modules::PanGenome::PanGenome;
use Modules::Phylogeny::PhylogenyFileCreator;
use Parallel::ForkManager;
use Tie::Log4perl;
use Log::Log4perl;
use File::Path 'make_path';
use Carp;
use Archive::Zip;
use Role::Tiny::With;

with 'Roles::CombineFilesIntoSingleFile';

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


=head3 _initialize

Initializes the logger.
Assigns all values to class variables.
Anything else that the _initialize function does.

=cut

sub _initialize{
	my($self)=shift;

  

    my $settingsObj = shift;
   
    #get and parse the configuration file for Panseq
	$self->settings($settingsObj);
	$self->_createDirectories();

	#logging
    $self->logger(Log::Log4perl->get_logger()); 
	$self->logger->info("Logger initialized in Modules::Setup::Panseq");  
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

=head3 settings

All Panseq settings from config file.

=cut

sub settings{
	my $self=shift;
	$self->{'_settings'}=shift // return $self->{'_settings'};
}

sub run{
	my $self=shift;

	#if mode is set to loci, launch loci_finder, rather than Panseq
	if($self->settings->runMode eq 'loci'){
		$self->_launchLociFinder();
	}
	else{
		$self->_launchPanseq();
	}

	$self->_cleanUp();
	$self->_createZipFile();
}

=head2 _launchLociFinder

Runs the loci_finder.pl script.

=cut

sub _launchLociFinder{
	my $self = shift;

	
}


=head2 _launchPanseq

Runs either pan-genome or novel region finding, based on the $self->settings->runMode option

=cut

sub _launchPanseq{
	my $self = shift;

	#get the query/reference files to use in the comparison
	my $files = Modules::Setup::PanseqFiles->new(
		'queryDirectory'=>$self->settings->queryDirectory,
		'referenceDirectory'=>$self->settings->referenceDirectory // undef
	);

	#get novel regions
	my $novelIterator = Modules::NovelRegion::NovelIterator->new(
		'queryFile'=>$files->singleQueryFile($self->settings->baseDirectory . 'singleQueryFile.fasta'),
		'referenceFile'=>$files->singleReferenceFile($self->settings->baseDirectory . 'singleReferenceFile.fasta'),
		'panGenomeFile'=>$self->settings->baseDirectory . 'panGenome.fasta',
		'novelRegionsFile'=>$self->settings->baseDirectory . 'novelRegions.fasta',
		'settings'=>$self->settings
	);
	$novelIterator->run();

	$self->logger->info("Panseq mode set as " . $self->settings->runMode);
	#perform pan-genomic analyses
	if(defined $self->settings->runMode && $self->settings->runMode eq 'pan'){
		my $panObj = $self->_performPanGenomeAnalyses($files,$novelIterator);
		$self->_createTreeFiles($panObj->panGenomeOutputFile,$panObj->coreSnpsOutputFile);
	}
}


=head3

Given the $self->settings->baseDirectory,
test to see if it exists.
If it does, stop program with a message.
If it does not, create directory, as well as directory/logs for logging output.

=cut

sub _createDirectories{
	my $self=shift;

	if(-d $self->settings->baseDirectory){
		print "\nThe directory you have specified for program output:\n\t" . $self->settings->baseDirectory . "\n"
			. 'already exists. Please specify a new baseDirectory' . "\n";
		exit(1);
	}
	else{
		#with File::Path
		make_path($self->settings->baseDirectory);
		make_path($self->settings->baseDirectory . 'logs/');
	}
}


=head3 _cleanUp

Removes any files left-over from the run.
Uses Roles::CombineFilesIntoSingleFile for _getFileNamesFromDirectory

=cut

sub _cleanUp{
	my $self=shift;
	
	my $fileNamesRef = $self->_getFileNamesFromDirectory($self->settings->baseDirectory);

	foreach my $file(@{$fileNamesRef}){
		if(
			($file =~ m/_db(\.n|temp)/) || 
			($file =~ m/_sequenceSplitterDB/) || 
			($file =~ m/singleQueryFile\.fasta/) ||
			($file =~ m/nucmer\.delta/) ||
			($file =~ m/Temp_\d\d\d\d\d/) ||
			($file =~ m/\.xml/)
		){
			unlink $file;
		}
	}
}

=head3 _createTreeFiles

Using the Modules::Phylogeny::PhylogenyFileCreator, create both a core-snp
and pan-genome +/- file for use in phylogenetic analyses.
If more than one processor available, make both at the same time.

=cut

sub _createTreeFiles{
	my $self = shift;
	my $panGenomeFile=shift;
	my $coreSnpsFile=shift;

	my $forker= Parallel::ForkManager->new($self->settings->numberOfCores);

	for my $num(1..2){
		$forker->start and next;
			if($num==1){
				my $treeMaker=Modules::Phylogeny::PhylogenyFileCreator->new(
					'inputFile'=>$coreSnpsFile,
					'outputFileType'=>'phylip',
					'outputFile'=>$self->settings->baseDirectory . 'core_snps.phylip'
				);
				$treeMaker->run();
			}
			elsif($num==2){
				my $treeMaker=Modules::Phylogeny::PhylogenyFileCreator->new(
					'inputFile'=>$panGenomeFile,
					'outputFileType'=>'phylip',
					'outputFile'=>$self->settings->baseDirectory . 'pan_genome.phylip',
					'conversionFile'=>$self->settings->baseDirectory . 'phylip_name_conversion.txt'
				);
				$treeMaker->run();
			}
			else{
				$self->logger->logconfess("num value is $num, should not exceed 2");
			}
		$forker->finish();
	}
	$forker->wait_all_children();
}

sub _performPanGenomeAnalyses{
	my $self=shift;
	my $files=shift;
	my $novelIterator=shift;

	#fragmentationSize defaults to 0 if no fragmentation is to be done

	my $segmenter = Modules::Fasta::SegmentMaker->new(
		'inputFile'=>$novelIterator->panGenomeFile,
		'outputFile'=>$self->settings->baseDirectory . 'pangenome_fragments.fasta',
		'segmentSize'=>$self->settings->fragmentationSize
	);

	if($self->settings->fragmentationSize > 0){
		$segmenter->segmentTheSequence;
	}

	my $dbCreator = Modules::Alignment::MakeBlastDB->new(
		'dbtype'=>'nucl',
		'blastDirectory'=>$self->settings->blastDirectory,
		'in'=>$files->singleQueryFile,		
		'out'=>$self->settings->baseDirectory . 'panseq_db',
		'title'=>'panseq_db',
		'logfile'=>'>>'.$self->settings->baseDirectory . 'logs/FormatBlastDB.pm.log'	
	);
	$dbCreator->run();

	 #'out' can be a directory or file. If numberOfSplits
	 #is greater than 1, that many cores will be used and the query
	 #will be split approximately evenly among the processes.
	 #'out' will be used as the base from which the actual output files
	 #are created
	my $blaster = Modules::Alignment::BlastRun->new(
		'query'=>$segmenter->outputFile,
		'blastDirectory'=>$self->settings->blastDirectory,
		'task'=>'blastn',
		'db'=>$self->settings->baseDirectory . $dbCreator->title,
		'outfmt'=>'5',
		'evalue'=>'0.00001',
		'word_size'=>20,
		'num_threads'=>1,
		'numberOfSplits'=>20,
		'out'=>$self->settings->baseDirectory  
	);
	$blaster->run();
	#do the pan-genome analysis

	my $panAnalyzer = Modules::PanGenome::PanGenome->new(
		'xmlFiles'=>$blaster->outputXMLfiles,
		'numberOfCores'=>$self->settings->numberOfCores,
		'resultArraySize'=>5,
		'percentIdentityCutoff'=>$self->settings->percentIdentityCutoff,
		'coreGenomeThreshold'=>$self->settings->coreGenomeThreshold,
		'outputDirectory'=>$self->settings->baseDirectory,
		'muscleExecutable'=>$self->settings->muscleExecutable,
		'accessoryType'=>$self->settings->accessoryType,
		'queryFile'=>$files->singleQueryFile
	);
	$panAnalyzer->run();
	return $panAnalyzer;
}


sub _createZipFile{
    my $self=shift;
    
    $self->logger->info("Creating zip file");
    #object to add all files the user receives
    my $zipper = Archive::Zip->new();
    
    #this is used in creating the batch file for the program run. called 'baseDirectory' within the Panseq scripts.
    my $outputDir =  $self->settings->baseDirectory;

    #get all file names from directory (to avoid dealing with logging, not using FileManipulation.pm)
    opendir( DIRECTORY, $outputDir ) or die "cannot open directory $outputDir $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;
    foreach my $fileName (@dir) {
            #dont add directories or . and .. files
            next if substr( $fileName, 0, 1 ) eq '.';
            next if ($fileName eq 'logs');                  

            #usage is: addFile( $fileName [, $newName ] ) from Archive::Zip manual
            $zipper->addFile($outputDir . $fileName);
    }
    my $status = $zipper->writeToFileNamed($outputDir . 'panseq_results.zip');
    
    if($status == 0){
            #everything ok
            return 1;
    }
    else{
            #something went wrong, return error
            return 0;
    }
}




#do the core/accessory analysis



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

1;
