#!/usr/bin/env perl

package Modules::Setup::CombineFilesIntoSingleFile;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../";
use IO::File;

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

    $self->logger->debug("Logger initialized in Modules::Setup::CombineFilesIntoSingleFile");  
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
    my $self=shift;
    $self->{'_logger'} = shift // return $self->{'_logger'};
}



=head3 _combineFilesIntoSingleFile

Takes all filenames from the arrayref passed to it and prints out to the filename given to it.
Both are required options and the program will die if either are absent.
Uses the ->getlines function of IO::File rather than printing line by line.
This should not be a memory concern, but will change if it turns out to be.
The append flag allows adding to a file, rather than erasing the contents and making a new output file.
If append is not specified, it defaults to off.

=cut

sub combineFilesIntoSingleFile{
	my $self=shift;
	my $filesToCombine = shift // $self->logger->logdie("_combineFilesIntoSingleFile requires files in an array ref");
    my $outputFile = shift // $self->logger->logdie("outputFile required in _combineFilesIntoSingleFile $!");
    my $append = shift // 0;
    my $firstLineAdjust = shift // 0;

    my $filemode;
    if($append == 1){
       $filemode = '>>';
    }
    else{
        $filemode = '>';
    }

    open(my $coutFH, $filemode, $outputFile) or $self->logger->fatal("Could not open $outputFile\n $!");
    $self->logger->debug("Combining files:");

    my $firstFile = 1;
	foreach my $file(@{$filesToCombine}){
        unless(-e $file && -s $file > 0){
            $self->logger->warn("Skipping $file as it has size of 0");
            next;
        }

        $self->logger->debug("Adding $file to $outputFile");

        open(my $cinFH, '<', $file) or $self->logger->fatal("Cannot open $file $!");

        if($firstLineAdjust && ($firstFile > 1)){
            $cinFH->getline();
        }    

        while(my $cline = $cinFH->getline()){
            $cline =~ s/\R//g;
            $coutFH->print($cline . "\n");
        }
        $cinFH->close();
        $firstFile++;
	}
    $coutFH->close();
    return $outputFile;
}

sub combineAndSanitizeFastaFiles{
    my $self=shift;

    my $filesToCombine = shift // $self->logger->logdie("_combineAndSanitizeFastaFiles requires files in an array ref.");
    my $outputFile = shift // $self->logger->logdie("outputFile required in _combineAndSanitizeFastaFiles $!");

    my $outFH = IO::File->new('>'. $outputFile) or $self->logger->logdie("Could not open $outputFile\n $!");

    $self->logger->debug("Combining files:");
    foreach my $file(@{$filesToCombine}){
        $self->logger->debug("$file");
        my $inFH = IO::File->new('<' . $file) or die "Cannot open $file\n $!";
        
        while(my $line = $inFH->getline){
            $line =~ s/\R//g;
            if($line =~ m/^>/){
                $line =~ s/[^\w\|\.\=\>]/_/g;                
            }
            $line .= "\n";
            $outFH->print($line);
        }

        $inFH->close();
    }
    $outFH->close();
    return $outputFile;
}

=head3 _getFileNamesFromDirectory

Opens the specified directory, excludes all filenames beginning with '.' and
returns the rest as an array ref.

=cut

sub getFileNamesFromDirectory{
    my $self=shift;
    my $directory = shift;

    opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;

    my @fileNames;
    foreach my $fileName(@dir){
        next if substr( $fileName, 0, 1 ) eq '.';
        push @fileNames, ( $directory . $fileName );
    }

    return \@fileNames;
}

1;
