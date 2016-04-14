#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use File::Copy;

my $genomeList = $ARGV[0] // die "no genome list file specified!\n";
my $htmlFile = $ARGV[1] // die "no HTML file specified!\n";


#make a copy of the original, read it in
#write out over the original, with the new data
#keep the temp file around for disasters

#with File::Copy
my $tempFile = $htmlFile . '_temp';
copy($htmlFile, $tempFile) or die "Could not create temp html file!\n";


#write out the new html file with new data
open(my $genomeFH, '<', $genomeList) or die "$!";
open(my $htmlFH, '>', $htmlFile) or die "$!";


my %genomeHash=();
while(my $line = $genomeFH->getline()){
    next if $. == 1;

    $line =~ s/\R//g;
    my @la = split('\t', $line);

    my $organismName = $la[7];
    my $strainName = $la[8];

    my $ftpName;
    if($la[19] =~ m/genomes\/all\/(.+)/){
        $ftpName= $1;
    }
    else{
        die "Could not get FTP location for $organismName\n";
    }

    $genomeHash{$ftpName}=$organismName . ', ' . $strainName;
}

open(my $tempFH, '<', $tempFile) or die "Could not open $tempFile $!\n";
my $switch = 0;
my $prefix = 'q_';
while(my $line = $tempFH->getline()){
    if($switch == 1){
        #skip all the existing list elements
        if($line =~ m/\<li\>/){
            next;
        }
        else{
            foreach my $k(sort {$genomeHash{$a} cmp $genomeHash{$b}} keys %genomeHash){
                $htmlFH->print('<li><input type="checkbox" name="' . $prefix . $k  . '">' . $genomeHash{$k} . '</li>' . "\n");
            }
            $switch=0;
        }
    }
    $htmlFH->print($line);

    if($line =~ m/\Flip the switch/){
        $switch = 1;
    }
}

