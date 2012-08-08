#!/usr/bin/perl

=head1 PangenomeGenerator.pm

groupSelection.pl - based on coreAccessory.pm code, this module simply generates a non-redudant pangenome

=head1 SYNOPSIS

    my $ca = PangenomeGenerator->new();
	$ca->PangenomeGenerator($configFile);
    END

=head1 DESCRIPTION

Uses a batchfile as input. To run see example batchfile for formatting.

=cut

package PangenomeGenerator;

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
use Carp;
use FileInteraction::FileManipulation;
use NovelRegion::NovelRegionFinder;
use FileInteraction::Fasta::SequenceRetriever;
use Pipeline::PanseqShared;
use FileInteraction::FlexiblePrinter;
use TreeBuilding::PhylogenyFileCreator;
use Data::Dumper;
use Logging::Logger;
use Object::Tiny::RW qw{
  _segmentCoreInput
  _novelSeedName
  _seedFileName
  _notSeedFileName
  _fragmentationSize
};

our @ISA = qw/PanseqShared FlexiblePrinter/;

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_coreAccessoryInitialize(@_);
	return $self;
}

#methods
sub _coreAccessoryInitialize {
	my ($self) = shift;

	#inheritance
	$self->_panseqSharedInitialize(@_);
	$self->_flexiblePrinterInitialize(@_);
}

#methods

sub _createNovelConfigFile {
	my ($self) = shift;

	my $fileName = $self->_baseDirectory . 'temp_novel_config.txt';
	my $outFH = IO::File->new( '>' . $fileName ) or die "$!";
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

sub _createSeedAndNotSeedFiles {
	my ($self) = shift;

	my $seed = $self->_getSeedFromQueryNameHash();

	#create database from queryFiles
	my $retriever = SequenceRetriever->new( $self->combinedQueryFile );
	$self->_seedFileName( $self->_baseDirectory . 'seedFile_' . $seed->name . '.fasta' );
	$self->_notSeedFileName( $self->_baseDirectory . 'notSeedFile' . '.fasta' );

	my $seedFH    = IO::File->new( '>' . $self->_seedFileName )    or die "Cannot open " . $self->_seedFileName . "$!";
	my $notSeedFH = IO::File->new( '>' . $self->_notSeedFileName ) or die 'Cannot open ' . $self->_notSeedFileName . "$!";

	#create reference file from seed
	foreach my $queryName ( keys %{ $self->queryNameObjectHash } ) {
		foreach my $fastaHeader ( @{ $self->queryNameObjectHash->{$queryName}->arrayOfHeaders } ) {
			if ( $self->queryNameObjectHash->{$queryName}->name eq $seed->name ) {
				print $seedFH '>' . $fastaHeader . "\n" . $retriever->extractRegion($fastaHeader) . "\n";
			}
			else {
				print $notSeedFH '>' . $fastaHeader . "\n" . $retriever->extractRegion($fastaHeader) . "\n";
			}
		}
	}
	$seedFH->close();
	$notSeedFH->close();
}

sub runPangenomeGenerator {
	my ($self) = shift;

	if (@_) {
		my $configFile = shift;

		#initialization
		$self->_validateCoreSettings( $self->getSettingsFromConfigurationFile($configFile) );
		$self->createDirectories();

		my $inputLociFile;
		$self->getQueryNamesAndCombineAllInputFiles();

		$self->logger->debug( "DEBUG:\tAfter combining all input files there are: " . ( scalar keys %{ $self->queryNameObjectHash } ) );
		$self->logger->info( "INFO:\t" . 'Determining non-redundant pan-genome' );
		$self->_createSeedAndNotSeedFiles();

		#get no_duplicates novel regions
		my $novelRegionFinder = NovelRegionFinder->new();

		#runNovelRegionFinder return file name of novelRegions file
		my $novelRegionFile = $novelRegionFinder->runNovelRegionFinder( $self->_createNovelConfigFile );

		#combine the novel regions with the seed file for a complete pan-genome
		my $combiner = FileManipulation->new();
		$inputLociFile = $self->_baseDirectory . 'panGenome.fasta';
		my $combinedFH = IO::File->new( '>' . $inputLociFile );
		$combiner->outputFilehandle($combinedFH);
		$combiner->vanillaCombineFiles( [ ( $self->_seedFileName, $novelRegionFile ) ] );
		$combinedFH->close();

		#segment the input if required
#		if ( $self->_segmentCoreInput eq 'yes' ) {
			$self->logger->info("INFO:\tSegmenting the input sequences");

			#create outputFH
			my $segmentedFile = $self->_baseDirectory . 'inputLoci_segments.fasta';
			my $segmentedPanFH = IO::File->new( '>' . $segmentedFile ) or die "$!";

			my $segmenter = SegmentMaker->new($segmentedPanFH);
			$segmenter->segmentTheSequence( $inputLociFile, $self->_fragmentationSize );
			$segmentedPanFH->close();
#		}
	}
}

sub _validateCoreSettings {
	my ($self) = shift;

	my $validator = $self->_validator;

	if (@_) {
		my $settingsHashRef = shift;

		foreach my $setting ( keys %{$settingsHashRef} ) {
			my $value = $settingsHashRef->{$setting};

			$self->_segmentCoreInput( $validator->yesOrNo($value) )  if $setting eq 'segmentCoreInput';
			$self->_fragmentationSize( $validator->isAnInt($value) ) if $setting eq 'fragmentationSize';
		}
	}
}

sub _getSeedFromQueryNameHash {
	my ($self) = shift;

	my $seed;
	foreach my $queryName ( keys %{ $self->queryNameObjectHash } ) {
		$seed = $queryName;
		my $numberOfFastaSequences = scalar( @{ $self->queryNameObjectHash->{$queryName}->arrayOfHeaders } );

		#if single fasta header, means a closed genome
		if ( $numberOfFastaSequences == 1 ) {
			last;
		}
	}
	$self->logger->info( "INFO:\t" . $self->queryNameObjectHash->{$seed}->name . ' selected as seed sequence.' );

	return $self->queryNameObjectHash->{$seed};
}

1;
