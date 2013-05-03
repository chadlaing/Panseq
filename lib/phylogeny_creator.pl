#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/";
use Modules::Phylogeny::PhylogenyFileCreator;
use Log::Log4perl qw/:easy/;
use Getopt::Long;

Log::Log4perl->easy_init($DEBUG);

#using Getopt::Long
my %options;
GetOptions(
	'inputFile=s'=>\$options{'inputFile'},
	'outputFile=s'=>\$options{'outputFile'},
	'conversionFile=s'=>\$options{'conversionFile'},
	'outputFormat:s'=>\$options{'outputFormat'},
	'help:s'=>\$options{'help'}
);

if(defined $options{'help'} || !defined $options{'inputFile'} || !defined $options{'outputFile'} || !defined $options{'conversionFile'}){
	printHelp();
}
else{
	my $phylo = Modules::Phylogeny::PhylogenyFileCreator->new(
		'inputFile'=>$options{'inputFile'},
		'outputFile'=>$options{'outputFile'},
		'conversionFile'=>$options{'conversionFile'},
		'outputFormat'=>$options{'outputFormat'}
	);
	$phylo->run();	
}

sub printHelp{
	print "\nphylogeny_creator.pl has the following options:\n\t",
	"--inputFile <file>\t:(required) specifies the tabular data to create a phylogeny from.\n\t",
	"--outputFile <file>\t:(required) specifies the location of the phylip formatted output file.\n\t",
	"--conversionFile <file>\t:(required) specifies the location of the name conversion file. Translates phylip number to original name.\n\t",
	"--outputFormat <phylip>\t:(optional) defaults to 'phylip'. Currently only phylip format file output is supported.\n\t",
	"--help\t\t\t:This help message.\n\n";
}

