#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';


my $inputDir = $ARGV[0];
my $fileNames = _getFileNamesFromDirectory($inputDir);


foreach my $file(@{$fileNames}){
    my $filename;

    if($file =~ m/(^.+)/){
        $filename = $1;
        $filename =~ s/\./_/g;
    }
    else{
        print STDERR "unknown file type: $filename";
        exit(1);
    }

    my $inFH = IO::File->new('<'. $inputDir . $file) or die "$!";
    my $outFH = IO::File->new('>' . $inputDir . $filename . '.checked') or die "$!";

    while(my $line = $inFH->getline){
        if($line =~ m/^>/){
            $line =~ s/>/>lcl\|$filename\|/;
        }
        $outFH->print($line);
    }
    $inFH->close();
    $outFH->close();

    #remove unchecked file
    unlink($inputDir . $file);
}



=head3 _getFileNamesFromDirectory

Opens the specified directory, excludes all filenames beginning with '.' and
returns the rest as an array ref.

=cut

sub _getFileNamesFromDirectory{
    my $directory = shift;

    opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;

    my @fileNames;
    foreach my $fileName(@dir){
        next if substr( $fileName, 0, 1 ) eq '.';
        next if -d $fileName;
        push @fileNames, ($fileName );
    }

    return \@fileNames;
}

