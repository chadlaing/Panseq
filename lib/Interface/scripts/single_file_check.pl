#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use FindBin;
use lib "$FindBin::Bin/../../";
use Modules::Fasta::SequenceName;


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
        $line =~ s/\R//g;
        my $sn = Modules::Fasta::SequenceName->new($line);


        #if the line is the same as the sequence name, then the header
        #is being used directly, which means multiple fasta sequences
        #will be treated as distinct genomes.
        #to overcome this, the filename is used as the lcl|| name.
        #this means people can upload files without modding the headers
        #if there is a single genome per file.
        if($line eq ('>' . $sn->name())){
            $line =~ s/>/>lcl\|$filename\|/;
        }
        $outFH->print($line . "\n");
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

