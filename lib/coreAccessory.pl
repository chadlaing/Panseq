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

my $configFile='/home/phac/workspace/Panseq_standalone/lib/core_config.txt';

#run
my $ca = CoreAccessory->new();
$ca->runCoreAccessory($configFile);