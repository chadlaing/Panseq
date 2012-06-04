#!/usr/bin/perl

package FileManipulation;

use strict;
use warnings;
use FindBin::libs;
use IO::File;
use Carp;
use FileInteraction::FlexiblePrinter;
use Log::Log4perl;
our @ISA = qw/FlexiblePrinter/;

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_fileManipulationInitialize(@_);
	return $self;
}

#class variables
sub maxLineSize {
	my $self = shift;
	$self->{'_fileManipulation_maxLineSize'} = shift // return $self->{'_fileManipulation_maxLineSize'};
}

sub logger{
	my $self = shift;
	$self->{'_fileManipulation_logger'} = shift // return $self->{'_fileManipulation_logger'};
}

sub fastaHeaders{
	my $self=shift;
	$self->{'_fileManipulation_fastaHeaders'} = shift // return $self->{'_fileManipulation_fastaHeaders'};
}

#methods
sub _fileManipulationInitialize {
	my $self = shift;

	#inheritance
	$self->_flexiblePrinterInitialize(@_);
	
	#logging
	$self->logger(Log::Log4perl->get_logger());

	#defaults
	$self->maxLineSize(10000);
	$self->fastaHeaders([]);
}

sub getFileNamesFromDirectory {
	my ($self) = shift;

	if (@_) {
		my $directory = shift;
		my @fileNames;
		opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
		my @dir = readdir DIRECTORY;
		closedir DIRECTORY;

		$self->logger->info("Getting file names from $directory");

		foreach my $fileName (@dir) {
			next if substr( $fileName, 0, 1 ) eq '.';

			#$self->logger->debug
			$self->logger->debug("Getting $fileName from $directory");

			#$self->logger->debug
			push @fileNames, ( $directory . $fileName );
		}
		return \@fileNames;
	}
	else {
		print STDERR "no directory specified!\n";
		exit(1);
	}
}

sub getFastaHeadersFromFile {
	my ($self) = shift;

	if (@_) {
		my $file    = shift;
		my $toClean = 0;
		if (@_) {
			$toClean = 1;
		}
		my @fastaHeaders;
		my $inFile = IO::File->new( '<' . $file ) or die "$!";

		while ( my $line = $inFile->getline ) {
			if ( $line =~ /^>(.+)/ ) {
				my $header = $1;
				$header = $self->cleanLine($header) if $toClean == 1;
				push @fastaHeaders, $header;
			}
		}
		$inFile->close;
		$self->fastaHeaders(\@fastaHeaders);
		return \@fastaHeaders;
	}
	else {
		print STDERR "nothing sent to getFastaHeadersFromFile\n";
		exit(1);
	}
}

sub combineFilesIntoSingleFile {
	my ($self) = shift;

	if ( scalar(@_) >= 1 ) {
		my $arrayRef = shift;

		my $toClean = 0;
		if (@_) {
			my $option = shift;
			$toClean = 1 if $option == 1;
		}

		foreach my $file ( @{$arrayRef} ) {

			#$self->logger->debug
			$self->logger->debug( "Adding $file to combined file: " . $self->outputFilehandle );

			#$self->logger->debug

			my $inFile = IO::File->new( '<' . $file ) or die "$! \n Cannot open $file \n";
			while ( my $line = $inFile->getline ) {
				if ($toClean) {
					$line = $self->cleanLine($line) if ( $line =~ /^>/ );
				}
				$self->printOut($line);
			}
			$self->printOut("\n");
			$inFile->close();
		}
	}
	else {
		print STDERR "wrong number of arguments sent to combineFilesIntoSingleFile!\n";
		exit(1);
	}
}

sub cleanLine {
	my ($self) = shift;

	my $line = shift;
	$line =~ s/[\s\'\"\.\:\+]/_/g;
	return ( $line . "\n" );
}

#setters
#Default is 10000 characters
sub setMaxLineSize {
	my $self = shift;
	my $size = shift;

	$self->maxLineSize($size);

	unless ( $self->maxLineSize
		&& ( $self->maxLineSize == int( $self->maxLineSize ) )
		&& ( $self->maxLineSize > 0 ) )
	{
		confess 'Max Line Size must be an integer greater than 0' . "\n";
	}
}

#getters
sub getMaxLineSize {
	my $self = shift;
	return $self->maxLineSize;
}

sub vanillaCombineFiles {
	my ($self) = shift;
	my $arrayRef = shift;
	my $destroy = shift // 0;

	$self->logger->debug("In vanillaCombineFiles");

	foreach my $file ( @{$arrayRef} ) {
		my $inFile = IO::File->new( '<' . $file ) or die "$! \n Cannot open $file \n";

		$self->logger->debug("combining $file");

		while ( my $line = $inFile->getline ) {
			$self->printOut($line);
		}
		$inFile->close();
		unlink $file if $destroy;
	}
}

sub combineDeltaFiles {
	my ($self) = shift;
	my $arrayRef = shift;
	my $referenceName = shift;
	my $queryName = shift;
	my $destroy = shift // 0;

	#header
	my $header = $referenceName . ' ' . $queryName . "\n" . 'NUCMER' . "\n";

	$self->logger->debug("In combineDeltaFiles");
	$self->printOut($header);

	foreach my $file ( @{$arrayRef} ) {
		my $inFile = IO::File->new( '<' . $file ) or die "$! \n Cannot open $file \n";

		$self->logger->debug("combining $file");

		#skip header lines (first 3 lines)
		print $inFile->getline();
		print $inFile->getline();
		

		while ( my $line = $inFile->getline ) {
			$self->printOut($line);
		}
		$inFile->close();
		unlink $file if $destroy;
	}
}
1;
