#!/usr/bin/perl

#Module allowing for the parallization of mummer given 1 query file with multiple fasta sequences
#Outwardly, the script will result in increased memory usage, processor usage (more than one)
#but will speed up processing time will yielding the same output as a single threaded Mummer
#be carful about number of cores used as memory usage effectively doubles at each core added
#Script intended for use with nucmer 3.22

package MummerParallel;

use strict;
use warnings;
use FindBin::libs;
use Carp;
use IO::File;
use File::Copy;
use Parallel::ForkManager;
use FileInteraction::Fasta::FastaFileSplitter;
use Mummer::MummerIO;
use Cwd 'getcwd';
use Logging::Logger;
our @ISA = qw/Logger/;

use Object::Tiny::RW qw{
  _b
  _c
  _d
  _g
  _l
  _p
  _queryFile
  _referenceFile
  _mummerDirectory
  _nooptimize
  _forkNumber
  _numberOfTempQueryFiles
  _arrayOfTempQueryFiles
  _arrayOfTempDeltaFiles
  deltaFile
  baseDirectory
};

#constructor; Usage: MummerParallel->new($self->_mummerDirectory, $self->_referenceFile, $self->_queryFile, $self->_forkNumber )
sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;

	$self->_mummerDirectory(shift);
	$self->_referenceFile(shift);
	$self->_queryFile(shift);
	$self->_forkNumber(shift);

	unless ( $self->_mummerDirectory
		&& $self->_queryFile
		&& $self->_referenceFile
		&& $self->_forkNumber )
	{
		confess
			'One or more required constructor parameters missing. Usage: MummerParallel->new($self->_mummerDirectory, $self->_referenceFile, $self->_queryFile, $self->_forkNumber )';
	}

	$self->_nooptimize(0);
	$self->_p('out');
	$self->_mummerParallelInitialize(@_);
	
	return $self;
}

#init
sub _mummerParallelInitialize{
	my $self=shift;
	
	$self->_loggerInitialize(@_);
	
	$self->_arrayOfTempQueryFiles([]); #anon
	$self->_arrayOfTempDeltaFiles([]); #anon	
	
	#to accomodate specifying a directory where MummerParallel should operate
	$self->baseDirectory(getcwd() . '/') unless defined $self->baseDirectory;
}

#setters
sub setB {
	my ($self) = shift;
	$self->_b( int(shift) );
	unless ( $self->_b && ( $self->_b >= 0 ) ) {
		confess 'b be an integer greater or equal to 0.';
	}
}

sub setC {
	my ($self) = shift;
	$self->_c( int(shift) );
	unless ( $self->_c && ( $self->_c >= 0 ) ) {
		confess 'c be an integer greater or equal to  0.';
	}
}

sub setD {
	my ($self) = shift;
	$self->_d(shift);
	unless ( $self->_d && ( $self->_d >= 0 ) ) {
		confess 'd be an float greater or equal to 0.';
	}
}

sub setG {
	my ($self) = shift;
	$self->_g( int(shift) );
	unless ( $self->_g && ( $self->_g >= 0 ) ) {
		confess 'g be an integer greater or equal to 0.';
	}
}

sub setL {
	my ($self) = shift;
	$self->_l( int(shift) );
	unless ( $self->_l && ( $self->_l >= 0 ) ) {
		confess 'l be an integer greater or equal to 0.';
	}
}

sub setP {
	my ($self) = shift;
	$self->_p(shift);
	unless ( $self->_p ) {
		confess 'setP must have input; usage: setP($name)';
	}
}

sub setNooptimizeOn {
	my ($self) = shift;
	$self->_nooptimize(1);
}

#getters
sub getB {
	my ($self) = shift;
	return $self->_b;
}

sub getC {
	my ($self) = shift;
	return $self->_c;
}

sub getD {
	my ($self) = shift;
	return $self->_d;
}

sub getG {
	my ($self) = shift;
	return $self->_g;
}

sub getL {
	my ($self) = shift;
	return $self->_l;
}

sub getP {
	my ($self) = shift;
	return $self->_p;
}

sub getQueryFile {
	my ($self) = shift;
	return $self->_queryFile;
}

sub getReferenceFile {
	my ($self) = shift;
	return $self->_referenceFile;
}

sub getMummerDirectory {
	my ($self) = shift;
	return $self->_mummerDirectory;
}

sub isNooptimizeOn {
	my ($self) = shift;
	return $self->_nooptimize;
}

#primary methods
sub run {
	my ($self) = shift;
	
	$self->logger->info(
		"INFO:\tStarting MummerParallel run\n\t" . 
		'Mummer Directory: ' . $self->_mummerDirectory . "\n\t" .
		'Query File: ' . $self->_queryFile . "\n\t" .
		'Reference File: ' . $self->_referenceFile . "\n\t" .
		'Number of Processes: ' . $self->_forkNumber	
	);
	
	#run settings
	my %runSettingsHash=(
		'referenceFile' => $self->_referenceFile,
		'runDirectory'  => (getcwd() . '/'),
		'b'             => $self->_b,
		'c'             => $self->_c,
		'd'             => $self->_d,
		'g'             => $self->_g,
		'l'             => $self->_l			
	);			
	$self->deltaFile($self->baseDirectory . $self->_p . '.delta');	
	
	if ( $self->_forkNumber > 1 ) {
		my $forkManager     = Parallel::ForkManager->new( $self->_forkNumber );
		
		#create the split fasta files
		my $ffSplitter = FastaFileSplitter->new();
		$ffSplitter->splitFastaFile(
			$self->_queryFile,
			$self->_forkNumber
		);
		
		foreach my $tempQueryFileName(@{$ffSplitter->arrayOfSplitFiles}){
			push @{$self->_arrayOfTempDeltaFiles}, ($tempQueryFileName . '.mummerDeltaTemp.delta');
			$forkManager->start and next();
			
			#create temp copy of reference file to avoid races
			my $tempRefFileName = $tempQueryFileName . '.mummerRefFastaTemp';
			copy($self->_referenceFile, $tempRefFileName);
			
			my $mummerRun = MummerIO->new( $self->_mummerDirectory );
			if ( $self->_nooptimize ) {
				$runSettingsHash{'nooptimize'}=$self->_nooptimize;
			}
			
			$runSettingsHash{'queryFile'}= $tempQueryFileName;
			$runSettingsHash{'referenceFile'}= $tempRefFileName;
			$runSettingsHash{'p'}= $tempQueryFileName . '.mummerDeltaTemp';
			
			$mummerRun->runMummer(%runSettingsHash);
			unlink( $tempQueryFileName);
			unlink($tempRefFileName);
			$forkManager->finish();
		}
		
		$forkManager->wait_all_children();		
		$self->_combineTempDeltaFiles();
	}
	else {
		my $mummerRun = MummerIO->new( $self->_mummerDirectory );
		if ( $self->_nooptimize ) {
			$runSettingsHash{'nooptimize'}=$self->_nooptimize;
			$runSettingsHash{'queryFile'}=$self->_queryFile;
			$runSettingsHash{'p'}=$self->baseDirectory. $self->_p;
		}
		$mummerRun->runMummer(%runSettingsHash);
	}
}

sub _combineTempDeltaFiles{
	my $self=shift;
	
	$self->logger->info("INFO:\tCombining temp delta files");
	
	my $combiner = FileManipulation->new();
	my $combinedFH = IO::File->new('>' . $self->deltaFile) or die 'Cannot open ' . $self->deltaFile . "$!\n";
	$combiner->outputFilehandle($combinedFH);
	$combiner->vanillaCombineFiles($self->_arrayOfTempDeltaFiles,1); #flag destroys file
	$combinedFH->close();
}

1;
