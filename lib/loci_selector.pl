#!/usr/bin/perl
use strict;
use warnings;
use FindBin::libs;

use LociSelector::LociSelector;

unless(scalar(@ARGV)>=2){
	print STDERR "incorrect number of arguments to LociSelector\n",
	"Correct usage is: \n\t 
	perl lociSelector.pl <input file name> <loci number, or 'best'> > <output file name>";
	exit(1);
}

my $fileName=shift;
my $lociNumber=shift;

my $finder = LociSelector->new();
$finder->getBestLoci(
	$fileName,
	$lociNumber
);
