#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/";
use Modules::LociSelector::LociSelector;
use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init($DEBUG);

unless(scalar(@ARGV)>=2){
        print STDERR "incorrect number of arguments to LociSelector\n",
        "Correct usage is: \n\t 
        perl lociSelector.pl <input file name> <loci number, or 'best'> > <output file name>";
        exit(1);
}

my $finder = Modules::LociSelector::LociSelector->new(
	inputFile=>$ARGV[0],
	lociNumber=>$ARGV[1],
	maximizePod=>$ARGV[2] // 0
);

$finder->run();

