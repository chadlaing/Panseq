#!/usr/bin/perl

package PanseqShared;

use strict;
use warnings;
use diagnostics;
use FindBin::libs;
use Carp;
use File::Copy;
use FileInteraction::LinePackage;
use Pipeline::Validator;
our @ISA = qw/Validator/;

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_panseqSharedInitialize(@_);
	return $self;
}

#class variables
sub queryDirectory {
	my $self = shift;
	$self->{'_PanseqShared_queryDirectory'} = shift // return $self->{'_PanseqShared_queryDirectory'};
}

sub referenceDirectory {
	my $self = shift;
	$self->{'_PanseqShared_referenceDirectory'} = shift // return $self->{'_PanseqShared_referenceDirectory'};
}

sub _baseDirectory {
	my $self = shift;
	$self->{'_PanseqShared_baseDirectory'} = shift // return $self->{'_PanseqShared_baseDirectory'};
}

sub _numberOfCores {
	my $self = shift;
	$self->{'_PanseqShared_numberOfCores'} = shift // return $self->{'_PanseqShared_numberOfCores'};
}

sub mummerDirectory {
	my $self = shift;
	$self->{'_PanseqShared_mummerDirectory'} = shift // return $self->{'_PanseqShared_mummerDirectory'};
}

sub minimumNovelRegionSize {
	my $self = shift;
	$self->{'_PanseqShared_minimumNovelRegionSize'} = shift // return $self->{'_PanseqShared_minimumNovelRegionSize'};
}

sub createGraphic {
	my $self = shift;
	$self->{'_PanseqShared_createGraphic'} = shift // return $self->{'_PanseqShared_createGraphic'};
}

sub combinedQueryFile {
	my $self = shift;
	$self->{'_PanseqShared_combinedQueryFile'} = shift // return $self->{'_PanseqShared_combinedQueryFile'};
}

sub queryNameObjectHash {
	my $self = shift;
	$self->{'_PanseqShared_queryNameObjectHash'} = shift // return $self->{'_PanseqShared_queryNameObjectHash'};
}

sub _novel_nucB {
	my $self = shift;
	$self->{'_PanseqShared_novel_nucB'} = shift // return $self->{'_PanseqShared_novel_nucB'};
}

sub _novel_nucC {
	my $self = shift;
	$self->{'_PanseqShared_novel_nucC'} = shift // return $self->{'_PanseqShared_novel_nucC'};
}

sub _novel_nucD {
	my $self = shift;
	$self->{'_PanseqShared_novel_nucD'} = shift // return $self->{'_PanseqShared_novel_nucD'};
}

sub _novel_nucG {
	my $self = shift;
	$self->{'_PanseqShared_novel_nucG'} = shift // return $self->{'_PanseqShared_novel_nucG'};
}

sub _novel_nucL {
	my $self = shift;
	$self->{'_PanseqShared_novel_nucL'} = shift // return $self->{'_PanseqShared_novel_nucL'};
}

sub _validator {
	my $self = shift;
	$self->{'_panseqShared_validator'} = shift // return $self->{'_panseqShared_validator'};
}

#methods
sub _panseqSharedInitialize {
	my $self = shift;
	$self->_validator(Validator->new() );
}

sub getSettingsFromConfigurationFile {
	my ($self) = shift;

	if (@_) {
		my $fileName = shift;
		my $inFile = IO::File->new( '<' . $fileName ) or die "\nCannot open $fileName\n $!";
		my %settingsHash;

		while ( my $line = $inFile->getline ) {
			next unless $line =~ m/\t/;
			my $la = LinePackage->new($line);

			my $setting = $la->lineArray->[0] // 'no_value';
			my $value   = $la->lineArray->[1] // 'no_value';
			$settingsHash{$setting} = $value;
		}
		$inFile->close();

		$self->validateUniversalSettings( \%settingsHash );
		return \%settingsHash;
	}
	else {
		confess "no file specified in getSettingsFromConfigurationFile\n";
	}
}

sub missingParam {
	my ($self) = shift;

	print STDERR "Missing parameter $_[0], which is required for a NovelRegionFinder run!\n";
	exit(1);
}

sub validateUniversalSettings {
	my ($self) = shift;

	if (@_) {
		my $settingsHashRef = shift;

		foreach my $setting ( keys %{$settingsHashRef} ) {
			my $value = $settingsHashRef->{$setting};

			#directories
			$self->queryDirectory( $self->_validator->isADirectory($value) )     if $setting eq 'queryDirectory';
			$self->_baseDirectory( $self->_validator->isADirectory($value) )     if $setting eq 'baseDirectory';
			$self->mummerDirectory( $self->_validator->isADirectory($value) )    if $setting eq 'mummerDirectory';
			$self->referenceDirectory( $self->_validator->isADirectory($value) ) if $setting eq 'referenceDirectory';

			#integers
			$self->_numberOfCores( $self->_validator->isAnInt($value) )         if $setting eq 'numberOfCores';
			$self->minimumNovelRegionSize( $self->_validator->isAnInt($value) ) if $setting eq 'minimumNovelRegionSize';

			#yes or no
			$self->createGraphic( $self->_validator->yesOrNo($value) ) if $setting eq 'createGraphic';
		}

		#defaults and requirements
		$self->missingParam('queryDirectory')         unless ( defined $self->queryDirectory );
		$self->missingParam('baseDirectory')          unless defined $self->_baseDirectory;
		$self->missingParam('mummerDirectory')        unless defined $self->mummerDirectory;
		$self->_numberOfCores(1)                      unless defined $self->_numberOfCores;
		$self->missingParam('minimumNovelRegionSize') unless defined $self->minimumNovelRegionSize;
		$self->createGraphic('no')                    unless defined $self->createGraphic;
	}
}

sub getQueryNamesAndCombineAllInputFiles {
	my ($self) = shift;

	if (1) {

		#get list of the query files
		#we need to create a queryFiles and referenceFiles file
		my $fileManipulator = FileManipulation->new();
		my $queryFileName   = $self->_baseDirectory . 'queryFilesCombined.fasta';
		my $queryOutputFH   = IO::File->new( '>' . $queryFileName ) or die "$!";
		$fileManipulator->outputFilehandle($queryOutputFH);    #set filehandle

		#store for posterity
		$self->combinedQueryFile($queryFileName);
		$fileManipulator->combineFilesIntoSingleFile( $fileManipulator->getFileNamesFromDirectory( $self->queryDirectory ), 1 );                                                    #second argument cleans file
		$queryOutputFH->close();

		if ( $self->can('novelRegionFinderMode') ) {
			my $referenceFileName = $self->_baseDirectory . 'referenceFilesCombined.fasta';

			$self->combinedReferenceFile($referenceFileName);    #save for posterity
			my $refFH = IO::File->new( '>' . $referenceFileName ) or die "$!";
			$fileManipulator->outputFilehandle($refFH);
			$fileManipulator->combineFilesIntoSingleFile( $fileManipulator->getFileNamesFromDirectory( $self->referenceDirectory ), 1 );
			$refFH->close();

			#for a unique novelRegionFinder run, need an all vs. all comparison
			my $outRefFH = IO::File->new( '>>' . $self->combinedReferenceFile ) or die "$!";
			$fileManipulator->outputFilehandle($outRefFH);
			$fileManipulator->vanillaCombineFiles( [ $self->combinedQueryFile ] );
			$outRefFH->close();

			#unique requires that all vs. all be performed
			if ( $self->novelRegionFinderMode eq 'unique' ) {

				#with File::Copy
				copy( $self->combinedReferenceFile, $self->combinedQueryFile ) or die "Could not complete File::Copy $!";
			}
		}

		$self->queryNameObjectHash( $self->getSequenceNamesAsHashRef( $fileManipulator->getFastaHeadersFromFile( $self->combinedQueryFile ) ) );
	}
	else {
		print STDERR "The universe has reversed polarity\n";
		exit(1);
	}
}

sub getSequenceNamesAsHashRef {
	my ($self) = shift;

	if (@_) {
		my $fastaHeadersArrayRef = shift;

		#get query names
		my %seqNames;
		foreach my $header ( @{$fastaHeadersArrayRef} ) {
			my $seqName = SequenceName->new($header);

			if ( defined $seqNames{ $seqName->name } ) {
				my $tempSeq = $seqNames{ $seqName->name };
				$tempSeq->addToArrayOfHeaders($header);
				$seqNames{ $seqName->name } = $tempSeq;
			}
			else {
				$seqNames{ $seqName->name } = $seqName;
			}
		}    #end foreach
		return \%seqNames;
	}    # end if
	else {
		print STDERR "nothing sent to getQuerySequenceNamesAsHashRef\n";
		exit(1);
	}
}

1;
