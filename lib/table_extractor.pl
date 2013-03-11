#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/";
use Modules::PanGenome::TableExtractor;
use Log::Log4perl qw/:easy/;
use Getopt::Long;

Log::Log4perl->easy_init($DEBUG);

#using Getopt::Long
my %options;
GetOptions(
	'inputFile=s'=>\$options{'inputFile'},
	'tableType=s'=>\$options{'tableType'},
	'minDiffChars:i'=>\$options{'minDiffChars'},
	'minPresent:i'=>\$options{'minPresent'},
	'maxMissing:i'=>\$options{'maxMissing'},
	'percentId:i'=>\$options{'percentId'},
	'help:s'=>\$options{'help'}
);

if(defined $options{'help'} || !defined $options{'inputFile'} || !defined $options{'tableType'}){
	printHelp();
}
else{
	my $tab = Modules::PanGenome::TableExtractor->new(
		'inputFile'=>$options{'inputFile'},
		'tableType'=>$options{'tableType'},
		'minimumCharacters'=>$options{'minDiffChars'} // 2,
		'maximumMissing'=>$options{'maxMissing'} // 0,
		'minimumPresent'=>$options{'minPresent'},
		'percentId'=>$options{'percentId'}
	);
	$tab->run();	
}

sub printHelp{
	print "\ntable_extractor.pl has the following options:\n\t",
	"-inputFile <file>\t\t:(required) specifies the tabular data to extract a subset from\n\t",
	"-tableType <core|pan>\t\t:(required) specifies either 'core' (SNP data) or 'pan' (%ID data)\n\t",
	#"-minimumCharacters <integer>\t:(optional) defaults to '2'. Specifies the minimum distinct character types required per line\n\t",
	"-maximumMissing <integer>\t:(optional) defaults to '0'. Specifies the maximum number of '-' characters allowed per line\n\t",
	"-minimumPresent <integer>\t:(optional) specifies the minimum number of non '-' characters required per line\n\t",
	"-percentId <integer>\t\t:(required for tableType 'pan') specifies the minimum %ID required for a locus to be considered 'present'\n\n";
}

