#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Modules::Fasta::SequenceRetriever;

my $sr = Modules::Fasta::SequenceRetriever->new();
$sr->printOut("Hi\n");
