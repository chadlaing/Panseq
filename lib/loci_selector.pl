#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/";

use Modules::LociSelector::LociSelector;

#_not_ready();


unless(scalar(@ARGV)>=2){
        print STDERR "incorrect number of arguments to LociSelector\n",
        "Correct usage is: \n\t 
        perl lociSelector.pl <input file name> <loci number, or 'best'> > <output file name>";
        exit(1);
}

my $fileName=shift;
my $lociNumber=shift;

my $finder = LociSelector::LociSelector->new();
$finder->getBestLoci(
        $fileName,
        $lociNumber
);

sub _not_ready{
	print STDERR "Still porting loci selector. It is currently non-functional\n";
	exit;
}