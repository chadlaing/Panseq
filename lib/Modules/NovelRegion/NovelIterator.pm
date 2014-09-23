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
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use File::Copy;
use Carp;
use Modules::Alignment::NucmerRun;
use Modules::NovelRegion::NovelRegionFinder;
use Parallel::ForkManager;
use Log::Log4perl;
use Digest::MD5;

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
    foreach my $key(sort keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::NovelRegion::NovelRegionFinder");
		}
	}

	#add the filenames for the hard-coded "last files" as the information cannot escape a fork
	$self->_lastNovelRegionsFile($self->settings->baseDirectory . 'lastNovelRegionsFile.fasta');
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

=head2 _lastNovelRegionsFile

Each iteration, the novel regions are output to this file.
The pangenome is built up by the iterative addition of the novelRegionsFile.

=cut

sub _lastNovelRegionsFile{
	my $self=shift;
	$self->{'_lastNovelRegions'}=shift // return $self->{'_lastNovelRegions'};
}


sub _retriever{
	my $self = shift;
	$self->{'_retriever'}=shift // return $self->{'_retriever'};
}


=head2 run

Launches the iterative novel region search.
For the pan-genome generation see _processRemainingThreeFiles description.

When running under 'unique' mode, both query and reference genomes as used for finding novel regions.
That way when running A vs. B, we end up with A(B) and B(A) from the single comparison. We also add
a reference genome to A(B) and B(A) so that the following comparisons will be against the totality of A + B.
For the reference [] indicates the fasta headers have been changed to lcl|_ignore|\w+(of previous header).
This instructs the no_duplicates method to ignore these as query, and not generate novel regions, which would
incorrectly generate non-unique sections. However, they are still used as reference, so that a unique query must
take into account all previous genomes.

Eg. 6 genomes (A,B,C,D,E,F), split among the cores:
A vs. B = A(B) + B(A) + [A]
C vs. D = C(D) + D(C) + [C]
E vs. F = E(F) + F(E) + [E]

A(B) + B(A) + [A] vs. C(D) + D(C) + [C]
	= A(BCD) + B(ACD) + C(DAB) + D(CAB) + [C(D) + D(C) + C]

If we have A + B vs. C, we will never get A vs. B.
To ensure this happens, run a self vs. self of the novel regions as a last step.

=cut

sub run{
	my $self=shift;

	my $numberOfGenomes = scalar(@{$self->settings->orderedGenomeNames});
	$self->logger->info("We have " . $numberOfGenomes . " genomes this run");
	if($numberOfGenomes ==0){
		$self->logger->fatal("Input error. Require at least 1 genome");
		exit(1);
	}

	#we need to run as many nucmer comparisons in parallel as possible
	#the total number of comparisons possible is <number of strains left> / 2
	#the total number of comparisons we can run at one time is in  $self->settings->numberOfCores
	#At the end, we return the non-redundant "pan-genome" from each core into a new list.
	#we iterate over this list until only a single file is left, the "pan-genome"

	$self->_retriever(Modules::Fasta::SequenceRetriever->new(
			'inputFile'=> $self->queryFile,
			'databaseFile'=>$self->settings->baseDirectory . 'queryfile_dbtemp'
		)
	);

	#create a file for each genome in the baseDirectory, for use in the nucmer comparisons
	my $allFastaFiles = $self->_createAllFastaFilesForNucmer($self->settings->orderedGenomeNames, $self->_retriever);
	my $numberOfRemainingFiles = scalar(@{$allFastaFiles});
	my $filesPerComparison=2;

	my $iterativeResultFile;
	while(scalar(@{$allFastaFiles}) > 1){
		$self->logger->info("Remaining files: " . scalar(@{$allFastaFiles}));

		$allFastaFiles = $self->_processRemainingFilesWithNucmer(
			$allFastaFiles,
			$filesPerComparison
		);

		$self->logger->debug("end of while:\n" . "@{$allFastaFiles}");
	}
	continue{
		if(scalar(@{$allFastaFiles})==1){
			$iterativeResultFile = $allFastaFiles->[0];
		}
	}	

	#double check
	unless(scalar(@{$allFastaFiles})==1){
		$self->logger->fatal("Failed to generate a single pan-genome file. Found " . scalar(@{$allFastaFiles}) . " files");
		exit(1);
	}

	my $finalFile = $iterativeResultFile;
	#need to run against the referenceDirectory files if they exist
	if($self->settings->novelRegionFinderMode eq 'unique'){
		$self->logger->info("Gathering unique novel regions.");

		my $finalReferenceFile;
		if(defined $self->referenceFile && -s $self->referenceFile > 0){		
			$self->logger->info("Comparing novel regions against reference file");
			
			my $combiner = Modules::Setup::CombineFilesIntoSingleFile->new();
			$finalReferenceFile= $combiner->combineFilesIntoSingleFile(
				[$self->_lastNovelRegionsFile, $self->referenceFile],
				$self->settings->baseDirectory . 'finalReferenceFile_unique.fasta'
			);
		}
		else{
			$self->logger->info("Using _lastNovelRegionsFile, as there is no reference file specified");
			$finalReferenceFile = $self->_lastNovelRegionsFile();
		}	

		my $coordsFile = $self->_processNucmerQueue($self->_lastNovelRegionsFile,$finalReferenceFile, $self->settings->baseDirectory . 'unique_regions');
		my $novelRegionsFile = $self->_printNovelRegionsFromQueue($coordsFile, $self->_lastNovelRegionsFile, ($iterativeResultFile . '_novelRegions'));
		$self->logger->info("Novel regions file $novelRegionsFile has size " . -s $novelRegionsFile);
		$finalFile = $novelRegionsFile;	
	}
	elsif($self->settings->novelRegionFinderMode eq 'no_duplicates'){
		if((defined $self->referenceFile) && (-s $self->referenceFile > 0)){
			$finalFile = $self->_processNucmerQueue($iterativeResultFile,$self->referenceFile,$iterativeResultFile . 'vs_refDir');
		}
	}
	else{
		$self->logger->fatal("Unknown novelRegionFinderMode: " . $self->novelRegionFinderMode);
		exit(1);
	}
	return $finalFile;
}


=head2 _processRemainingFilesWithNucmer

Takes in a list of fasta files that need to be compared among each other for novel regions.
Given the computed number of filesPerComparison, cycle through the files and
perform one round of nucmer comparisons. Return the generated "pan-genomes" from each set of comparisons,
which will be fed back into this sub until only one file is returned, which is the non-redundant pan-genome
for the original files.

We need to ensure the pan-genome has only one seed genome. For example.
6 genomes as query: (A,B,C,D,E,F).
In _processRemainingFilesWithNucmer we get:
A vs B = A(B) + B
C vs D = C(D) + D
E vs F = E(F) + F

If we were to allow one of these files to be the reference, and the rest query we would get the following:
A(B) + B + C(D) + D vs E(F) + F = A(BEF) + B(EF) + C(DEF) + D(EF)
Obviously, there is redundancy and un-compared genomes.
We need the following:

A(B) + B vs C(D) + D = A(BCD) + B(CD) + C(D) + D
A(BCD) + B(CD) + C(D) + D vs E(F) + F = A(BCDEF) + B(CDEF) + C(DEF) + D(EF) + E(F) + F

Where all have been compared, and there is no redundancy in un-compared regions.

8 genomes:
A vs B = A(B) + B
C vs D = C(D) + D
E vs F = E(F) + F
G vs H = G(H) + H

A(B) + B vs C(D) + D = A(BCD) + B(CD) + C(D) + Din
E(F) + F vs G(H) + H = E(FGH) + F(GH) + G(H) + H

A(BCD) + B(CD) + C(D) + D vs E(FGH) + F(GH) + G(H) + H
  = A(BCDEFGH) + B(CDEFGH) + C(DEFGH) + D(EFGH) + E(FGH) + F(GH) + G(H) + H 

=cut


sub _processRemainingFilesWithNucmer{
	my $self=shift;
	my $allFastaFiles = shift;
	my $filesPerComparison = shift;	

	my $forker = Parallel::ForkManager->new($self->settings->numberOfCores);

	my @filesToRun;	
	my $counter=1;
	my @outputFileNames;	

	my $fastaFile;
	my $numFilesMinusOne = scalar(@{$allFastaFiles})-1;
	foreach my $i (0 .. $numFilesMinusOne){
		$fastaFile = $allFastaFiles->[$i];
		push @filesToRun, $fastaFile;	
		unless(@filesToRun == $filesPerComparison){
			next;
		}
		
		#get digest of filename
		my $digester = Digest::MD5->new();
		my ($referenceFile, $queryFile) = @filesToRun;
		my $rf = $digester->add($referenceFile)->hexdigest;
		my $qf = $digester->add($queryFile)->hexdigest;

		my $newFileName = $self->settings->baseDirectory . $qf . '_' . $rf;
		push @outputFileNames, $newFileName;

		$forker->start and next;								
			my $coordsFile = $self->_processNucmerQueue($queryFile,$referenceFile, $newFileName);

			my $namesFile;
			if($self->settings->novelRegionFinderMode eq 'unique'){
				
				my $combiner = Modules::Setup::CombineFilesIntoSingleFile->new();			
				$namesFile = $combiner->combineFilesIntoSingleFile(
					[$queryFile,$referenceFile],
					$self->settings->baseDirectory . 'uniqueNovelRegions.fasta'
				);					
			}
			else{
				$namesFile = $queryFile;
			}

			my $novelRegionsFile = $self->_printNovelRegionsFromQueue($coordsFile, $namesFile, ($newFileName . '_NR'));	
			
			#special case for mode unique
			$self->logger->debug("combinedFile inputs:\nnrf: $novelRegionsFile\nrf : $referenceFile\nnfn: $newFileName");
			my $combinedFile = $self->_combineNovelRegionsAndReferenceFile($novelRegionsFile,$referenceFile,$newFileName);	
			$self->logger->info("Combined file: $combinedFile");			

			#remove temp files
			unlink $coordsFile;
			unlink $newFileName . '_filtered.delta';
			unlink $newFileName . '_query_dbtemp.index';				

			#keep the last novelRegionsFile under new name
			#with File::Copy
			copy($novelRegionsFile,$self->_lastNovelRegionsFile) or die "$!";
			#unlink $queryFile;
			#unlink $referenceFile;					
		$forker->finish;		
	}
	continue{
		$counter++;

		if(scalar(@filesToRun) == $filesPerComparison){
			@filesToRun=();
		}
		else{
			if($i == $numFilesMinusOne){
				$self->logger->debug("Adding $fastaFile to outputFileNames");
				push @outputFileNames, $fastaFile;
				$self->logger->debug("@outputFileNames");
			}
		}
	}
	$forker->wait_all_children;
	return \@outputFileNames;
}

=head2 _combineNovelRegionsAndReferenceFile

Takes the novel regions file and the reference file and combines them.
This serves as the "pan-genome" from the particular set of comparisons that were run.
When the entire process iterates down to a single file, we have the pan-genome for all the initial genomes.
When the mode is set to 'unique', we need the pan-genome to be built, but want to prevent the added
reference file from being used to generate new novel regions; thus we replace the header with lcl|_ignore|\w+

=cut

sub _combineNovelRegionsAndReferenceFile{
	my $self = shift;
	my $novelRegionsFile = shift;
	my $referenceFile = shift;
	my $outputFile = shift;

	my $refFileToAdd;
	if($self->settings->novelRegionFinderMode eq 'unique'){
		$refFileToAdd = $self->_setRefFileToIgnore($referenceFile);
	}
	else{
		$refFileToAdd = $referenceFile;
	}

	my $combiner = Modules::Setup::CombineFilesIntoSingleFile->new();
	$combiner->combineFilesIntoSingleFile(
		[$novelRegionsFile, $refFileToAdd],
		$outputFile
	);
	return $outputFile;
}


=head2 _setRefFileToIgnore

Need to ensure the reference file is not used for generating new novel regions
under mode "unique". By setting the lcl|_ignore| flag, when used as a query, the novel
region finder will ignore it. When used as a reference, will still show the query matches.

=cut

sub _setRefFileToIgnore{
	my $self = shift;
	my $refFileOriginal = shift;

	my $fileName = $refFileOriginal . '_ignore';
	my $refInFH = IO::File->new('<' . $refFileOriginal) or die "$!";
	my $refOutFH = IO::File->new('>' . $fileName) or die "$!";
	while(my $line = $refInFH->getline){
		if($line =~ m/^>/){
			$line =~ s/\W/_/g;
			$line = '>lcl|_ignore|' . $line . "\n";
		}
		$refOutFH->print($line);
	}
	$refInFH->close();
	$refOutFH->close();
	return $fileName;
}

=head2 _printNovelRegionsFromQueue

Takes in the coords file, the queryFile, and an outputFile.
The output file name is also taken in and returned from the sub.

=cut

sub _printNovelRegionsFromQueue{
	my $self = shift;
	my $coordsFile = shift;
	my $queryFile = shift;
	my $outputFile = shift;

	my $nrf = Modules::NovelRegion::NovelRegionFinder->new(
		'settings'=>$self->settings,
		'coordsFile'=>$coordsFile,
		'queryFile'=>$queryFile,
	);
	$nrf->findNovelRegions();

	$self->logger->info("File for db construction: $queryFile");
	my $databaseFile = $queryFile . '_dbtemp';

	#need to include reference file if $self->type eq 'unique'
	#can't just use the $self->retriever, as it does not include the sequence (4..5000) in the headers
	#we could do some math to account for those however, but not right now
	#right now we just want it to function.

	my $retriever = Modules::Fasta::SequenceRetriever->new(
		'inputFile'=> $queryFile,
		'databaseFile'=>$databaseFile
	);
	$nrf->printNovelRegions($retriever, $outputFile);

	unlink $databaseFile;
	return $outputFile;
}

=head2 _processNucmerQueue

Takes in a query and reference file, runs nucmer and returns the name of the
generated coords file. Includes all files in the "reference directory" if any are present.

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
		'percentIdentityCutoff'=>$self->settings->percentIdentityCutoff,
		'logFile'=> $outputFile . 'mummerlog.txt'
	);	
	$self->logger->info("Nucmer query file: $queryFile");
	$self->logger->debug("Nucmer ref file: $referenceFile");
	$nucmer->run();	
	unlink $outputFile . '.delta';
	unlink $outputFile . '.index';
	return $nucmer->coordsFile;	
}



=head2 _getNumberOfComparisonsToRun

Based on the number of cores and number of remaining files to process,
determines how many files should be processed per concurrent mummer instance, given that we
need at least 2 files per comparison. 

=cut

sub _getNumberOfFilesPerComparison{
	my $self=shift;
	my $numberOfFiles = shift;

	my $numFilesPerComp = int( $numberOfFiles / $self->settings->numberOfCores );

	#Must have at least 2 files per comparison
	if($numFilesPerComp < 2){
		$numFilesPerComp = 2;
	}

	return ($numFilesPerComp);
}



=head2 _createAllFastaFilesForNucmer

This extracts all fasta sequences associated with a particular genome to their own individual
file on the system. These files are then used by nucmer for parallel processing.
Also creates a single fasta file of the reference sequences if they exist.

=cut

sub _createAllFastaFilesForNucmer{
	my $self = shift;
	my $genomeNames=shift;
	my $retriever = shift;
	
	my @outputFiles;
	foreach my $genome(@{$genomeNames}){
		my $outputFile = $self->_getCleanFileName($genome);
		
		$self->_createFileFromGenomeName(
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
	return($time . int(rand(1000000)) . int(rand(1000000)));
}


=head2 _getBestSeed

This looks through all the sequences available from the _multiFastaSN object, and selects either a closed genome
or a genome that has the fewest contigs. Returns the name of the best seed.

=cut

sub _getBestSeed{
	my $self=shift;
	my $retriever = shift;

	my %numberOfContigs;

	foreach my $contig (keys %{$self->settings->_genomeNameFromContig}){
		if(defined $numberOfContigs{$self->settings->_genomeNameFromContig->{$contig}}){
			$numberOfContigs{$self->settings->_genomeNameFromContig->{$contig}} += 1;
		}
		else{
			$numberOfContigs{$self->settings->_genomeNameFromContig->{$contig}} =1;
		}
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
		'genome'=>$seedName,
		'retriever'=>$retriever,
		'outputFile'=>$self->panGenomeFile
	);

'retriever' is a Modules::Fasta::SequenceRetriever object that is already loaded with the
multiple fasta query file. 

'outputFile' is the file to write the output to.

'genome' is the name of the sequence that will be extracted from retriever.
Iterates over all fasta headers in the object.

=cut

sub _createFileFromGenomeName{
	my $self=shift;

	my %params = @_;

	#check the params
	unless(defined $params{'genome'}
		&& defined $params{'retriever'}
		&& defined $params{'outputFile'}
	){
		$self->logger->logconfess('_createFileFromGenomeName requires the following:
			"genome"=><name>,
			"retriever"=><obj>,
			"outputFile"=><name>
		');
	}

	my $outFH = IO::File->new('>' . $params{'outputFile'}) or die "Cannot open " . $params{'outputFile'} . "$!\n";
	foreach my $header(@{$self->settings->contigNamesFromGenome->{$params{'genome'}}}){
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
