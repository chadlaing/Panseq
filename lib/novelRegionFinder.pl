#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use FindBin::libs;
use NovelRegion::NovelRegionFinder;

#usage: novelRegionFinder.pl <config file>


my $configFile;
if(defined $ARGV[0]){
	$configFile=$ARGV[0];
}
else{
	print STDERR "Please specify a configuration file!\n";
	exit(1);
}


#run
my $nrf = NovelRegionFinder->new();
$nrf->runNovelRegionFinder($configFile);
