#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Net::FTP;
use Archive::Extract;

my $localNcbiDirectory = $ARGV[0] // die "local NCBI directory required";

#sets up the parameters for ncbi ftp connection
my $host = 'ftp.ncbi.nlm.nih.gov';

#new FTP directory structure
my $directory = '/genomes/genbank/bacteria/';

#within this directory is organism_name/all_assembly_versions/name/name.fna.gz
#eg. ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/Escherichia_coli/all_assembly_versions/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz

#constructs the connection
my $ftp = Net::FTP->new($host, Debug => 1,Passive => 1, Timeout => 1000) or die "Cannot connect to genbank: $@";
#log in as anonymous, use email as password
$ftp->login("anonymous",'chadlaing@inoutbox.com') or die "Cannot login ", $ftp->message;
$ftp->cwd($directory) or die "Cannot change working directory ", $ftp->message;

#get list of all bacterial species names
my @speciesList = $ftp->ls();

#go through each species in turn and check if the file has already been downloaded
#if it has not, download
my $counter = 0;
foreach my $species(@speciesList){
    print STDERR "species: $species\n";
    my $speciesDirectory = $directory . $species . '/all_assembly_versions/';
    $ftp->cwd($speciesDirectory) or die "Could not change directory " . $ftp->message();
    my @genomesList = $ftp->ls();

    #each genome has its own directory
    #get the .fna.gz file
    foreach my $genome(@genomesList){
        print STDERR "genome: $genome\n";

#        $counter++;
#        if($counter >= 5){
#            print STDERR "Found 5, finishing\n";
#            exit(1);
#        }

        my $genomeDirectory = $speciesDirectory . $genome . '/';
        $ftp->cwd($genomeDirectory) or die $ftp->message();

        #get the .fna.gz file
        my $genomeFile;
        my @fileList = $ftp->ls();

        foreach my $file(@fileList){
            if($file =~ m/_genomic\.fna\.gz/){
                $genomeFile = $file;
                last;
            }
        }

        unless(defined $genomeFile){
            print STDERR "Could not find a .gz file for $genome in $genomeDirectory\n";
            next;
        }

        my $localFile = $localNcbiDirectory . $genomeFile;
        my $decompressedLocalFile;

        if($localFile =~ m/(.+)\.gz/){
            $decompressedLocalFile = $1;
        }
        else{
            print STDERR "Could not find .gz file for $genome\n";
            next;
        }

        #check if the filename currently exists in the local directory
        if(-e $decompressedLocalFile){
            print STDERR "Current genome $genomeFile exists, skipping\n";
        }
        else{
            #ascii mode messes up the data, must use binary
            $ftp->binary();
            $ftp->get($genomeFile, $localFile) or die ("download of $genomeFile failed\n" and next);
            $ftp->ascii();

            #extract the compressed file
            my $extracter = Archive::Extract->new('archive'=>$localFile);
            $extracter->extract('to'=>$localNcbiDirectory) or die ("Could not extract $localFile\n" and next);
            unlink $localFile;
        }


    }
}

