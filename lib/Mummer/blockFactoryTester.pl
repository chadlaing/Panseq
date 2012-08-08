#!/usr/local/bin/perl

use FindBin::libs;
use warnings;
use strict;
use diagnostics;
use DeltaBlockFactory;
use Benchmark::Timer;
use IO::File;

#timer for benchmarking
my $t;
my $fact;

#my $output = new IO::File('>'.'MummerMsaSpeed.txt');

$t = Benchmark::Timer->new( skip => 0 );
for ( 0 .. 0 ) {
	$t->start(' Factory');
	my $fh = new IO::File(
		'<' . '/home/justin/Panseq_dev/Panseq2/lib/Mummer/input/out.delta' );
	$fact = new DeltaBlockFactory($fh);
	my $currentBlock = $fact->nextDeltaBlock();
	print $currentBlock->header;
	while ($currentBlock) {
		print $currentBlock->refName . ' '
		  . $currentBlock->queryName . ' '
		  . $currentBlock->refLength . ' '
		  . $currentBlock->queryLength
		  . "\n";
		  print $currentBlock->getPercentID() ."\n";
		foreach ( @{$currentBlock->gapArray} ) {
			print $_ . "\n";
		}
		print "0\n";
		$currentBlock = $fact->nextDeltaBlock();
	}

	$t->stop;
	undef($fact);
}
print( $t->report );

