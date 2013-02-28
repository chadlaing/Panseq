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


=head3 _initialize

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

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

=head3 settings

An object containing all of the settings from the user.

=cut

sub settings{
	my $self=shift;
	$self->{'_settings'}=shift // return $self->{'_settings'};
}

=head3 queryFile

The file name containing all the query sequences.

=cut

sub queryFile{
	my $self=shift;
	$self->{'_queryFile'}=shift // return $self->{'_queryFile'};
}

=head3 referenceFile

The file name containing all the reference sequences.

=cut

sub referenceFile{
	my $self=shift;
	$self->{'_referenceFile'}=shift // return $self->{'_referenceFile'};
}


=head3 panGenomeFile

The file that contains the built-up pan-genome.
The _getSeedName method looks for a closed sequence or one with the fewest contigs, and this will be output
to the panGenomeFile to use as a seed

=cut

sub panGenomeFile{
	my $self=shift;
	$self->{'_seedFile'}=shift // return $self->{'_seedFile'};
}

=head3 novelRegionsFile

Each iteration, the novel regions are output to this file.
The pangenome is built up by the iterative addition of the novelRegionsFile.

=cut

sub novelRegionsFile{
	my $self=shift;
	$self->{'_novelRegionsFile'}=shift // return $self->{'_novelRegionsFile'};
}


=head3 run

Launches the iterative novel region search.

=cut

sub run{
	my $self=shift;

	#generate a hash of all query sequences
	my $multiFastaSN = Modules::Fasta::MultiFastaSequenceName->new(
		'fileName'=>$self->queryFile
	);

	#ensure we can retrieve the fasta sequences from the query file
	my $retriever = Modules::Fasta::SequenceRetriever->new(
		'inputFile'=> $self->queryFile,
		'outputFile'=> $self->settings->baseDirectory . 'novelRegions.fasta'
	);

	#get all names of genomes from the queryFile
	my @genomeNames = map {$_} (keys %{$multiFastaSN->sequenceNameHash});

	#create the seed file and get the seed name
	#we only need to do this if there is no reference file provided
	#otherwise, the reference file can be the start of the 'panGenome'
	my $seedName='';
	unless(defined $self->referenceFile){
		$seedName = $self->_getBestSeed($multiFastaSN, $retriever);
		$self->_createFileFromGenomeName(
			'mfsn'=>$multiFastaSN,
			'genome'=>$seedName,
			'retriever'=>$retriever,
			'outputFile'=>$self->panGenomeFile
		);	
	}else{
		#specify a '1' flag for append in the Roles::CombineFilesIntoSingleFile method
		#order: filesToCombne (array ref), outputFile, appendFlag
		$self->_combineFilesIntoSingleFile(
			[$self->referenceFile],
			$self->panGenomeFile
		);
	}

	#loop through all names, running mummer and expanding the pangenome / all novelRegions file
	foreach my $genome(@genomeNames){
		if($genome eq $seedName){
			$self->logger->debug("Skipping seed $seedName in NovelIterator");
			next;
		}

		my $queryFile = $self->_getCleanFileName($genome);

		#create the query file for mummer
		$self->_createFileFromGenomeName(
			'mfsn'=>$multiFastaSN,
			'genome'=>$genome,
			'retriever'=>$retriever,
			'outputFile'=>$queryFile
		);	

		#run mummer
		my $nucmer = Modules::Alignment::NucmerRun->new(
			'b'=>200,
			'c'=>65,
			'd'=>0.12,
			'g'=>90,
			'l'=>20,
			'coordsFile'=>$self->settings->baseDirectory . 'nucmer.coords',
			'p'=>$self->settings->baseDirectory . 'nucmer',
			'mummerDirectory'=>$self->settings->mummerDirectory,
			'queryFile'=>$queryFile,
			'referenceFile'=>$self->_getReferenceFileForNucmer(),
			'logFile'=>$self->settings->baseDirectory . 'mummerlog.txt'
		);	
		$self->logger->info("Nucmer query file: $queryFile");
		$nucmer->run();	

		#gather novel regions from coords file
		my $nrf = Modules::NovelRegion::NovelRegionFinder->new(
			'mode'=>$self->settings->novelRegionFinderMode,
			'coordsFile'=>$nucmer->coordsFile,
			'queryFile'=>$nucmer->queryFile,
			'minimumNovelRegionSize'=>$self->settings->minimumNovelRegionSize
		);
		$nrf->findNovelRegions();

		my $novelOutputFile = $self->settings->baseDirectory . $self->_getTempName();
		$nrf->printNovelRegions($retriever, $novelOutputFile);

		#build up the pan-genome
		#specify a '1' flag for append in the Roles::CombineFilesIntoSingleFile method
		#order: filesToCombne (array ref), outputFile, appendFlag
		$self->_combineFilesIntoSingleFile(
			[$novelOutputFile],
			$self->panGenomeFile,
			1
		);

		#also create a novelRegionsFile
		#this is useful when a reference file is specified
		#we only want the novelty from each query sequence, not the reference file as well
		$self->_combineFilesIntoSingleFile(
			[$novelOutputFile],
			$self->novelRegionsFile,
			1
		);
		#cleanup of files
		unlink $novelOutputFile;
		unlink $nucmer->queryFile;
		unlink $nucmer->coordsFile;
	}#end of foreach
}


=head3 _getTempName

Returns a filename based on the system time, for use as temp files.

=cut

sub _getTempName{
	my $self=shift;

	my $time = localtime();
	$time =~ s/\W//g;
	return($time . 'temp');
}


=head3 _getBestSeed

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


=head3 _createFileFromGenomeName

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


=head3 _getCleanFileName

Removes all pipes (|) and quotes ('',"") and returns the correct baseDirectory prefix for the query file.

=cut
sub _getCleanFileName{
	my $self=shift;
	my $name = shift;

	$name =~ s/[\|\'\"]/_/g;
	return $self->settings->baseDirectory . $name;
}

=head3 _getReferenceFileForNucmer

We need to include the growing no_duplicates (panGenomeFile) as well as any specified reference files
in the referenceFile. Creates a blank panGenomeFile if one does not already exist.

=cut

sub _getReferenceFileForNucmer{
	my $self = shift;

	#either a pan-genome file needs to exist with the seed, or a reference file with sequences needs to exist
	if(!defined $self->referenceFile && !(-e $self->panGenomeFile)){
		$self->logger->logconfess("Either a reference file or pan-genome file required in Modules::NovelIterator::_getReferenceFileForNucmer");
	}

	#create an empty file for panGenomeFile if it does not yet exist
	unless(-e $self->panGenomeFile){
		my $panFH = IO::File->new('>' . $self->panGenomeFile) or $self->logger->logdie("Could not open " . $self->panGenomeFile);
		$panFH->close();
	}


	if (defined $self->referenceFile){
		$self->logger->info('Combining ' . $self->referenceFile . ' and ' . $self->panGenomeFile . ' as single reference file for nucmer');

		return $self->_combineFilesIntoSingleFile(
			[$self->referenceFile,$self->panGenomeFile],
			$self->settings->baseDirectory . 'nucmerReferenceFile.fasta'
		);
	}
	else{
		return $self->panGenomeFile;
	}
}



1;
