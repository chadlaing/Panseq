#!/usr/bin/env perl

=pod

=head1 NAME

Modules::NovelRegion::NovelIterator - Iterates through the query genomes, building up a no_duplicates file
(pan-genome). TODO:: Iterative common to all.

=head1 SYNOPSIS

	use Modules::NovelRegion::NovelIterator;
	
	my $novelIterator = Modules::NovelRegion::NovelIterator->new(
		'queryFile'=>$files->singleQueryFile($settings->baseDirectory . 'singleQueryFile.fasta'),
		'panGenomeFile'=>$settings->baseDirectory . 'panGenome.fasta',
		'settings'=>$settings
	);
	$novelIterator->run();

=head1 DESCRIPTION

This module requires a Modules::Setup::Settings object to be passed, as well as a queryFile containing all the 
query sequences in multi-fasta format.
The panGenomeFile is the 'no_duplicates' output file that is iteratively created.

=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing(chadlaing gmail com)

=head2 Methods

=cut

package Modules::NovelRegion::NovelIterator;

#includes
use strict;
use warnings;
use diagnostics;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Carp;
use Modules::Alignment::NucmerRun;
use Modules::NovelRegion::NovelRegionFinder;
use Modules::Fasta::MultiFastaSequenceName;
use Parallel::ForkManager;
use Log::Log4perl;
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


=head2 _initialize

Initializes the logger.
Assigns all values to class variables.
Anything else that the _initialize function does.

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 

    $self->logger->debug("Logger initialized in Modules::NovelRegion::NovelIterator");  

    my %params = @_;

    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::NovelRegion::NovelRegionFinder");
		}
	}	
}

=head2 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

=head2 settings

An object containing all of the settings from the user.

=cut

sub settings{
	my $self=shift;
	$self->{'_settings'}=shift // return $self->{'_settings'};
}

=head2 queryFile

The file name containing all the query sequences.

=cut

sub queryFile{
	my $self=shift;
	$self->{'_queryFile'}=shift // return $self->{'_queryFile'};
}

=head2 referenceFile

The file name containing all the reference sequences.

=cut

sub referenceFile{
	my $self=shift;
	$self->{'_referenceFile'}=shift // return $self->{'_referenceFile'};
}


=head2 panGenomeFile

The file that contains the built-up pan-genome.
The _getSeedName method looks for a closed sequence or one with the fewest contigs, and this will be output
to the panGenomeFile to use as a seed

=cut

sub panGenomeFile{
	my $self=shift;
	$self->{'_seedFile'}=shift // return $self->{'_seedFile'};
}

=head2 novelRegionsFile

Each iteration, the novel regions are output to this file.
The pangenome is built up by the iterative addition of the novelRegionsFile.

=cut

sub novelRegionsFile{
	my $self=shift;
	$self->{'_novelRegionsFile'}=shift // return $self->{'_novelRegionsFile'};
}


=head2 run

Launches the iterative novel region search.

=cut

sub run{
	my $self=shift;

	#generate a hash of all query sequences
	my $multiFastaSN = Modules::Fasta::MultiFastaSequenceName->new(
		'fileName'=>$self->queryFile
	);

	#get all names of genomes from the queryFile
	my @genomeNames = map {$_} (keys %{$multiFastaSN->sequenceNameHash});

	my $numberOfGenomes = scalar(@genomeNames);
	$self->logger->info("We have " . $numberOfGenomes . " genomes this run");
	if($numberOfGenomes ==0){
		$self->logger->fatal("Input error. Require at least 1 genome");
		exit(1);
	}

	#we need to run as many nucmer comparisons in parallel as possible
	#the total number of comparisons possible is <number of strains left> / 2
	#the total number of comparisons we can run at one time is in  $self->settings->numberOfCores
	#If there are more comparisons than cores, we add additional strains to each comparison group 
	#so that all strains are distributed across all cores.
	#At the end, we return the non-redundant "pan-genome" from each core into a new list.
	#we iterate over this list in the same manner until 3 or fewer "pan-genomes" remain, 
	#at which point a single-core comparison of the remaining strains returns the final pan-genome.

	my $retriever = Modules::Fasta::SequenceRetriever->new(
		'inputFile'=> $self->queryFile,
		'databaseFile'=>$self->settings->baseDirectory . 'queryfile_dbtemp'
	);

	#create a file for each genome in the baseDirectory, for use in the nucmer comparisons
	my $allFastaFiles = $self->_createAllFastaFilesForNucmer(\@genomeNames, $multiFastaSN, $retriever);
	my $numberOfRemainingFiles = scalar(@{$allFastaFiles});

	while($numberOfRemainingFiles > 1){
		$self->logger->info("Remaining files: $numberOfRemainingFiles");
		my ($filesPerComparison,$remainder) = $self->_getNumberOfFilesPerComparison($numberOfRemainingFiles);

		$self->logger->info("Files per: $filesPerComparison, remainder: $remainder");

		$allFastaFiles = $self->_processRemainingFilesWithNucmer(
			$allFastaFiles,
			$filesPerComparison,
			$remainder,
			$retriever
		);
	}
	continue{
		$numberOfRemainingFiles = scalar(@{$allFastaFiles});
	}
}


=head2 _processRemainingFilesWithNucmer

Takes in a list of fasta files that need to be compared among each other for novel regions.
Given the computed number of filesPerComparison and the remainder, cycle through the files and
perform one round of nucmer comparisons. Return the generated "pan-genomes" from each set of comparisons,
which will be fed back into this sub until only one file is returned, which is the non-redundant pan-genome
for the original files.

=cut


sub _processRemainingFilesWithNucmer{
	my $self=shift;
	my $allFastaFiles = shift;
	my $filesPerComparison = shift;
	my $remainder = shift;
	my $retriever = shift;

	$self->logger->info("In _processRemainingFilesWithNucmer. filesPerComparison: $filesPerComparison");

	my $forker = Parallel::ForkManager->new($self->settings->numberOfCores);

	my @filesToRun;
	my @filesFromNucmer;

	my $counter=1;
	foreach my $fastaFile (@{$allFastaFiles}){

		
		if((!defined $filesToRun[0]) || (scalar(@filesToRun) < $filesPerComparison)){
			$self->logger->info("Pushing to array, size of:" . scalar(@filesToRun));
			push @filesToRun, $fastaFile;
		}
		else{
			$self->logger->info("In the else");
			if($remainder > 0){
				push @filesToRun, $fastaFile;
				$remainder--;
			}
			my $newFileName = $self->settings->baseDirectory . 'nucmerTempFile' . $counter . $self->_getTempName;
			push @filesFromNucmer, $newFileName;

			$forker->start and next;
				my ($queryFile, $referenceFile) = $self->_getQueryReferenceFileFromList(\@filesToRun);
				my $coordsFile = $self->_processNucmerQueue($queryFile,$referenceFile, $newFileName);
				my $novelRegionsFile = $self->_printNovelRegionsFromQueue($coordsFile, $queryFile, ($newFileName . 'novelRegions'), $retriever);	
				
				#combines them into the $newFileName specified above, which is added to the @filesFromNucmer array
				#this array is returned and used to feed back into this processing sub
				$self->_combineNovelRegionsAndReferenceFile($novelRegionsFile,$referenceFile,$newFileName);	

				#remove temp files
				unlink $queryFile;
				unlink $referenceFile;
				unlink $coordsFile;
				unlink $novelRegionsFile;
			$forker->finish;
		}
	}
	continue{
		$counter++;
		if(scalar(@filesToRun) > $filesPerComparison){
			$self->logger->info("Resetting filesToRun from size of " . scalar(@filesToRun));
			@filesToRun=();
		}
	}

	$forker->wait_all_children;
	return \@filesFromNucmer;
}


=head2 _combineNovelRegionsAndReferenceFile

Takes the novel regions file and the reference file and combines them.
This serves as the "pan-genome" from the particular set of comparisons that were run.
When the entire process iterates down to a single file, we have the pan-genome for all the initial genomes.

=cut

sub _combineNovelRegionsAndReferenceFile{
	my $self = shift;
	my $novelRegionsFile = shift;
	my $referenceFile = shift;
	my $outputFile = shift;

	#with Roles::CombineFilesIntoSingleFile
	$self->_combineFilesIntoSingleFile(
		[$novelRegionsFile, $referenceFile],
		$outputFile
	);
	return $outputFile;
}


=head2 _printNovelRegionsFromQueue

Takes in the coords file, the queryFile, and a retriever object.
The retriever has the database of all initial query sequences and is passed around
to avoid having to regenerate the database.
The output file name is also taken in and returned from the sub.

=cut

sub _printNovelRegionsFromQueue{
	my $self = shift;
	my $coordsFile = shift;
	my $queryFile = shift;
	my $outputFile = shift;
	my $retriever = shift;

	my $nrf = Modules::NovelRegion::NovelRegionFinder->new(
		'mode'=>$self->settings->novelRegionFinderMode,
		'coordsFile'=>$coordsFile,
		'queryFile'=>$queryFile,
		'minimumNovelRegionSize'=>$self->settings->minimumNovelRegionSize
	);
	$nrf->findNovelRegions();
	$nrf->printNovelRegions($retriever, $outputFile);
	return $outputFile;
}


=head2 _getQueryReferenceFileFromList

Takes in a list of file names.
Makes the first item of the list the reference file.
Combines the remaining files into a single "queryFile".
Returns the name of this queryFile.
We combine query files, as nucmer "streams" this file past the reference file,
so larger query files are better than larger reference files.
The individual files combined into the single query file are deleted.

=cut

sub _getQueryReferenceFileFromList{
	my $self = shift;
	my $listOfFiles = shift;

	my $ref = shift @{$listOfFiles};

	#with Roles::CombineFilesIntoSingleFile
	my $queryFileName = $self->settings->baseDirectory . 'query' . $self->_getTempName;
	$self->_combineFilesIntoSingleFile(
		$listOfFiles,
		$queryFileName
	); 

	#remove the non-combined query files
	foreach my $file(@{$listOfFiles}){
		unlink $file;
	}
	return ($queryFileName,$ref);
}


=head2 _processNucmerQueue

Takes in a query and reference file, runs nucmer and returns the name of the
generated coords file.

=cut

sub _processNucmerQueue{
	my $self = shift;
	my $queryFile = shift;
	my $referenceFile = shift;
	my $outputFile = shift;

	#run mummer
	my $nucmer = Modules::Alignment::NucmerRun->new(
		'b'=>200,
		'c'=>65,
		'd'=>0.12,
		'g'=>90,
		'l'=>20,
		'coordsFile'=>$outputFile . '.coords',
		'p'=>$outputFile,
		'mummerDirectory'=>$self->settings->mummerDirectory,
		'queryFile'=>$queryFile,
		'referenceFile'=>$referenceFile,
		'logFile'=> $outputFile . 'mummerlog.txt'
	);	
	$self->logger->info("Nucmer query file: $queryFile");
	$nucmer->run();	
	return $nucmer->coordsFile;	
}



=head2 _getNumberOfComparisonsToRun

Based on the number of cores and number of remaining files to process,
determines how many files should be processed per concurrent mummer instance, given that we
need at least 2 files per comparison. Returns the number of files that are the
"remainder" and must be added to a group.

=cut

sub _getNumberOfFilesPerComparison{
	my $self=shift;
	my $numberOfFiles = shift;

	my $numFilesPerComp = int( $numberOfFiles / $self->settings->numberOfCores );

	#Must have at least 2 files per comparison
	if($numFilesPerComp < 2){
		$numFilesPerComp = 2;
	}

	my $remainder = $numberOfFiles % $numFilesPerComp;

	return ($numFilesPerComp,$remainder);
}



=head2 _createAllFastaFilesForNucmer

This extracts all fasta sequences associated with a particular genome to their own individual
file on the system. These files are then used by nucmer for parallel processing.
Also creates a single fasta file of the reference sequences if they exist.

=cut

sub _createAllFastaFilesForNucmer{
	my $self = shift;
	my $genomeNames=shift;
	my $mfsn=shift;
	my $retriever = shift;
	
	my @outputFiles;
	foreach my $genome(@{$genomeNames}){
		my $outputFile = $self->_getCleanFileName($genome);
		
		$self->_createFileFromGenomeName(
			'mfsn'=>$mfsn,
			'genome'=>$genome,
			'retriever'=>$retriever,
			'outputFile'=>$outputFile
		);

		push @outputFiles, $outputFile;
	}

	return \@outputFiles;
}



=head2 _getTempName

Returns a filename based on the system time, for use as temp files.

=cut

sub _getTempName{
	my $self=shift;

	my $time = localtime();
	$time =~ s/\W//g;
	return($time . 'temp');
}


=head2 _getBestSeed

This looks through all the sequences available from the _multiFastaSN object, and selects either a closed genome
or a genome that has the fewest contigs. Returns the name of the best seed.

=cut

sub _getBestSeed{
	my $self=shift;
	my $mfsn = shift;
	my $retriever = shift;

	my %numberOfContigs;

	foreach my $genome (keys %{$mfsn->sequenceNameHash}){
		$numberOfContigs{(scalar(@{$mfsn->sequenceNameHash->{$genome}->arrayOfHeaders}))}=$genome;
	}

	my $smallestGenome;
	foreach my $key(sort keys %numberOfContigs){
		$smallestGenome=$numberOfContigs{$key};
		last;
	}
	
	$self->logger->info("Best seed: $smallestGenome");
	return $smallestGenome;
}


=head2 _createFileFromGenomeName

Extracts all fasta sequences associated with a genome name from the queryFile containing all
query sequences in multi-fasta format.

Takes in named parameters as follows:
	->createFileFromGenomeName(
		'mfsn'=>$multiFastaSN,
		'genome'=>$seedName,
		'retriever'=>$retriever,
		'outputFile'=>$self->panGenomeFile
	);

'retriever' is a Modules::Fasta::SequenceRetriever object that is already loaded with the
multiple fasta query file. 

'mfsn' is a Modules::Fasta::MultiFastaSequenceName object, already loaded with the SequenceNames
of the input query file.

'outputFile' is the file to write the output to.

'genome' is the name of the sequence that will be extracted from retriever.
Iterates over all fasta headers in the object.

=cut

sub _createFileFromGenomeName{
	my $self=shift;

	my %params = @_;

	#check the params
	unless(defined $params{'mfsn'}
		&& defined $params{'genome'}
		&& defined $params{'retriever'}
		&& defined $params{'outputFile'}
	){
		$self->logger->logconfess('_createFileFromGenomeName requires the following:
			"mfsn"=><obj>,
			"genome"=><name>,
			"retriever"=><obj>,
			"outputFile"=><name>
		');
	}

	my $outFH = IO::File->new('>' . $params{'outputFile'}) or die "Cannot open " . $params{'outputFile'} . "$!\n";
	foreach my $header(@{$params{'mfsn'}->sequenceNameHash->{$params{'genome'}}->arrayOfHeaders}){
		$outFH->print('>' . $header . "\n" . $params{'retriever'}->extractRegion($header) . "\n");
	}
	$outFH->close();

	$self->logger->info("Creating " . $params{'outputFile'} . " from " . $ params{'genome'});
}


=head2 _getCleanFileName

Removes all pipes (|) and quotes ('',"") and returns the correct baseDirectory prefix for the query file.

=cut
sub _getCleanFileName{
	my $self=shift;
	my $name = shift;

	$name =~ s/[\|\'\"]/_/g;
	return $self->settings->baseDirectory . $name;
}


1;
