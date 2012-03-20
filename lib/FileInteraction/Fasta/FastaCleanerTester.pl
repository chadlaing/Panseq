#!/usr/local/bin/perl

use warnings;
use strict;
use diagnostics;
use FastaCleaner;
use Benchmark::Timer;
use IO::File;

#timer for benchmarking
my $t;
my $cleaner;

#my $output = new IO::File('>'.'MummerMsaSpeed.txt');

$t = Benchmark::Timer->new( skip => 0 );
for ( 0 .. 0 ) {
	$t->start(' MummerMSA');
	$cleaner =
	  new FastaCleaner( '/home/justin/Desktop/TestGenomes/senterica_draft_no_plasmids_100bp_limit.fasta','test.txt', 10000 );
	$cleaner->cleanup();

	$t->stop;
	undef($cleaner);
}
print( $t->report );

