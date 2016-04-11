#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use File::Copy;

my $dir = $ARGV[0];
my $htmlFile = $ARGV[1] // die "no HTML file specified!\n";

my $fileNamesRef = _getFileNamesFromDirectory($dir);

my %genomeList=();

foreach my $file(@{$fileNamesRef}){
    my $fullFile = $dir . $file;
    open(my $inFH, '<', $fullFile) or die "Cannot open $fullFile\n$!";

    while(my $line = $inFH->getline()){
        if($line =~ m/^>\S+\s+(.+?)(,|whole|assembly|genome|contig|complete|genomic|scaffold)/){
            my $name = $1;
            $genomeList{$name} = $file;
        }
        else{
            warn "Unable to parse $line\n";
        }
        last;
    }
    $inFH->close();
}

#make a copy of the original, read it in
#write out over the original, with the new data
#keep the temp file around for disasters

#with File::Copy
my $tempFile = $htmlFile . '_temp';
copy($htmlFile, $tempFile) or die "Could not create temp html file!\n";


#write out the new html file with new data
open(my $tempFH, '<', $tempFile) or die "$!";
open(my $htmlFH, '>', $htmlFile) or die "$!";

my $switch = 0;
while(my $line = $tempFH->getline()){
    if($switch == 1){
        #skip all the existing list elements
        if($line =~ m/\<li\>/){
            next;
        }
        else{
            foreach my $k(sort keys %genomeList){
                $htmlFH->print('<li><input type="checkbox" name="' . $genomeList{$k}  . '">' . $k . '</li>' . "\n");
            }
            $switch=0;
        }
    }
    $htmlFH->print($line);

    if($line =~ m/\Flip the switch/){
       $switch = 1;
    }
}



sub _getFileNamesFromDirectory{
    my $directory = shift;

    opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;

    my @fileNames;
    foreach my $fileName(@dir){
        next if substr( $fileName, 0, 1 ) eq '.';
        push @fileNames, $fileName;
    }

    return \@fileNames;
}
