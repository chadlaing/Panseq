#!/usr/bin/perl

use strict;
use warnings;
use FindBin::libs;
use Mummer::DeltaBlockFactory;
use IO::File;

my $FH = IO::File->new('<' . '/home/phac/workspace/Panseq_dev/Panseq2/lib/Testing/UnitTests/novelRegion/output_test/out.delta') or die "$!";

my $db=DeltaBlockFactory->new($FH);


while(my $block = $db->nextDeltaBlock()){
	print $block->refName . "\n",
		$block->refStart . "\n",
		$block->refEnd ."\n";
	
}

# DNA:
# Sequence 1:
#			ATGCGTACGACCCGGGTTT
#			TACGCATGCTGGGCCCAAA
#
# Sequence 2:
#			AAACCCGGGTCGTACGCAT
#			TTTGGGCCCAGCATGCGTA
