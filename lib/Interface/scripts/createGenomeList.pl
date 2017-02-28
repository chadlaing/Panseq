#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use File::Copy;

my $genomeList = $ARGV[0] // die "no genome list file specified!\n";
my $type = $ARGV[1] // die "no type of 'pan' or 'novel' specified!\n";



#we will make a _full and _reduced version to facilitate page loading
#_full will contain the entire set of genbank data (~15MB HTML file)
#_reduced version will contain every Nth item

#use the _base.html versions of novel and pan, to create the reduced and full
#version of each page
#the base version has no list elements, but is otherwise identical
#and includes the following text

my $REDUCED_TEXT_NOVEL = q{
    <p>
        IMPORTANT: To speed page loading, this view contains only 1% of the total genomes available for comparison. To load the complete set of genomes, please click <a href="/panseq/page/novel_full.html">Full Genome Set</a>
    </p>
};

my $REDUCED_TEXT_PAN = q{
    <p>
        IMPORTANT: To speed page loading, this view contains only 1% of the total genomes available for comparison. To load the complete set of genomes, please click <a href="/panseq/page/pan_full.html">Full Genome Set</a>
    </p>
};


#write out the new html file with new data
open(my $genomeFH, '<', $genomeList) or die "$!";

my %fullHash;
while(my $line = $genomeFH->getline()){
    next if $. == 1 || $. == 2;

    $line =~ s/\R//g;
    my @la = split('\t', $line);

    #19 is organism path
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
    if($la[19] =~ m/all\/(.+)$/){
        $ftpName= $1;
    }
    else{
        die "Could not get FTP location for $organismName from $la[19]\n";
    }
    $fullHash{$ftpName}=$fullName;
}

#go through base file
my $reducedFile;
my $fullFile;
my $baseFile;

if($type eq 'pan'){
    $reducedFile = 'pan.html';
    $fullFile = 'pan_full.html';
    $baseFile = 'pan_base.html';
}
elsif($type eq 'novel'){
    $reducedFile = 'novel.html';
    $fullFile = 'novel_full.html';
    $baseFile = 'novel_base.html';
}
else{
    die "Invalid type! Please specify novel or pan\n";
}

open(my $fullFH, '>', $fullFile) or die "$!";
open(my $reducedFH, '>', $reducedFile) or die "$!";

my $prefix;

open(my $baseFH, '<', $baseFile) or die "Couldn't open $baseFile $!\n";

while(my $line = $baseFH->getline()){
    if($line =~ m/--Reduced Only--/){
        if($type eq 'novel'){
            $reducedFH->print($REDUCED_TEXT_NOVEL);
        }
        elsif($type eq 'pan'){
            $reducedFH->print($REDUCED_TEXT_PAN);
        }
        else{
            die "Unknown type. Please choose novel or pan\n";
            exit(1);
        }
    }

    # <li><input type="checkbox" name="q_GCA_000160075.2_ASM16007v2">Abiotrophia defectiva ATCC 49176 </li>
    if($line =~m/--Query--/ || $line =~ m/--Reference--/){
        if($line =~ m/--Query--/){
            $prefix = 'q_';
        }
        else{
            $prefix = 'r_';
        }

        my $count=0;
        foreach my $k(sort {$fullHash{$a} cmp $fullHash{$b}} keys %fullHash){
            my $listLine = '<li><input type="checkbox" name="'. $prefix . $k . '">' . $fullHash{$k} . '</li>';

            #only print every 100 items to the reduced file
            if(($count % 100) == 0){
                $reducedFH->print($listLine . "\n");
            }
            $fullFH->print($listLine . "\n");

            $count++;
        }
    }
    else{
        $reducedFH->print($line);
        $fullFH->print($line);
    }
}
print "Finished\n";


