#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use File::Copy;

my $genomeList = $ARGV[0] // die "no genome list file specified!\n";
my $fullFile = $ARGV[1] // die "no full HTML file specified!\n";
my $reducedFile = $ARGV[2] // die "no reduced HTML file specified!";



#we will make a _full and _reduced version to facilitate page loading
#_full will contain the entire set of genbank data (~15MB HTML file)
#_reduced version will contain every Nth item
#make a copy of the original, read it in
#write out over the original, with the new data
#keep the temp file around for disasters

#with File::Copy
my $fullTempFile = $fullFile . '_temp';
my $reducedTempFile = $reducedFile . '_temp';
copy($fullFile, $fullTempFile) or die "Could not create temp full html file!\n";
copy($reducedFile, $reducedTempFile) or die "Could not create temp reduced html file\n";

#write out the new html file with new data
open(my $genomeFH, '<', $genomeList) or die "$!";
open(my $fullFH, '>', $fullFile) or die "$!";
open(my $reducedFH, '>', $reducedFile) or die "$!";


my %fullHash;
while(my $line = $genomeFH->getline()){
    next if $. == 1;

    $line =~ s/\R//g;
    my @la = split('\t', $line);

    my $organismName = $la[7];
    my $strainName = $la[8];
    my $isolateName = $la[9];
    my $asmName = $la[15];

    my $fullName = $organismName;
    if($strainName ne ""){
        $strainName =~ s/strain=//;

        #prevent duplicates like Campylobacter jejuni, Campylobacter jejuni,
        if($strainName ne $fullName){
            $fullName .= ', ' . $strainName;
        }
    }

    if($isolateName ne ""){
        $fullName .= ', ' . $isolateName;
    }

    if($asmName ne ""){
        $fullName .= ', ' . $asmName;
    }

    my $ftpName;
    if($la[19] =~ m/genomes\/all\/(.+)/){
        $ftpName= $1;
    }
    else{
        die "Could not get FTP location for $organismName\n";
    }
    $fullHash{$ftpName}=$fullName;
}

open(my $tempFH, '<', $reducedTempFile) or die "Could not open $reducedTempFile $!\n";
my $switch = 0;
my $reducedOnly = 0;
my $prefix = 'q_';
while(my $line = $tempFH->getline()){
    if($switch == 1){
        #skip all the existing list elements
        if($line =~ m/\<li\>/){
            next;
        }
        else{
            my $counter =0;
            foreach my $k(sort {$fullHash{$a} cmp $fullHash{$b}} keys %fullHash){
                if($counter % 100 == 0){
                    $reducedFH->print('<li><input type="checkbox" name="' . $prefix . $k  . '">' . $fullHash{$k} . '</li>' . "\n");
                }

                $fullFH->print('<li><input type="checkbox" name="' . $prefix . $k  . '">' . $fullHash{$k} . '</li>' . "\n");
                $counter++;
            }
            $switch=0;
        }
    }

    if($line =~ m/\<\!--Reduced Only--\>/){
        $reducedOnly++;
    }

    if($reducedOnly > 0){
        $reducedFH->print($line);
    }
    else{
        $fullFH->print($line);
        $reducedFH->print($line);
    }

    if($reducedOnly == 2){
        $reducedOnly=0;
    }


    if($line =~ m/\Flip the switch/){
        $switch = 1;
    }
}

