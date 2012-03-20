#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use FindBin::libs;
use MSA::BlastBased::CoreAccessory;
use Carp;

=head1
Takes in a single directory name and configuration file as input
Outputs all files to the directory specified in the configuration file.
=cut

my $configFile;
if(defined $ARGV[0]){
	$configFile=$ARGV[0];
}
else{
	print STDERR "Please specify a configuration file!\n";
	exit(1);
}


#run
my $ca = CoreAccessory->new();
$ca->runCoreAccessory($configFile);