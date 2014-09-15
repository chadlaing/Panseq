#!/usr/bin/env perl

package Roles::CombineFilesIntoSingleFile;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../";
use Role::Tiny;

=head3 _combineFilesIntoSingleFile

Takes all filenames from the arrayref passed to it and prints out to the filename given to it.
Both are required options and the program will die if either are absent.
Uses the ->getlines function of IO::File rather than printing line by line.
This should not be a memory concern, but will change if it turns out to be.
The append flag allows adding to a file, rather than erasing the contents and making a new output file.
If append is not specified, it defaults to off.

=cut

sub _combineFilesIntoSingleFile{
	my $self=shift;
	my $filesToCombine = shift // $self->logger->logdie("_combineFilesIntoSingleFile requires files in an array ref");
    my $outputFile = shift // $self->logger->logdie("outputFile required in _combineFilesIntoSingleFile $!");
    my $append = shift // 0;
    my $firstLineAdjust = shift // 0;

    my $outFH;

    if($append == 1){
        $outFH = IO::File->new('>>'. $outputFile) or $self->logger->logdie("Could not open $outputFile\n $!");
        $self->logger->debug("Append mode for $outputFile");
    }
    else{
        $outFH = IO::File->new('>'. $outputFile) or $self->logger->logdie("Could not open $outputFile\n $!");
    }
     
    $self->logger->debug("Combining files:");

    my $firstFile = 1;
	foreach my $file(@{$filesToCombine}){
        unless(-s $file > 0){
            $self->logger->info("Skipping $file as it has size of 0");
            next;
        }

        $self->logger->info("Adding $file to $outputFile");
        my $inFH = IO::File->new('<' . $file) or $self->logger->logdie("Cannot open $file $!");

        if($firstLineAdjust && $firstFile > 1){
            $self->logger->info("Discarding first line");
            $inFH->getline();
        }        
        $outFH->print($inFH->getlines);
        $inFH->close();

        $firstFile++;
	}
    $outFH->close();
    return $outputFile;
}

sub _combineAndSanitizeFastaFiles{
    my $self=shift;

    my $filesToCombine = shift // $self->logger->logdie("_combineAndSanitizeFastaFiles requires files in an array ref.");
    my $outputFile = shift // $self->logger->logdie("outputFile required in _combineAndSanitizeFastaFiles $!");

    my $outFH = IO::File->new('>'. $outputFile) or $self->logger->logdie("Could not open $outputFile\n $!");

    $self->logger->info("Combining files:");
    foreach my $file(@{$filesToCombine}){
        $self->logger->info("$file");
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

sub _getFileNamesFromDirectory{
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
