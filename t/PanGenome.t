#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Test::More 'no_plan';
use Test::Pretty;
use Log::Log4perl qw/get_logger/;
use Modules::PanGenome::PanGenome;

Log::Log4perl->easy_init();

my $pg = Modules::PanGenome::PanGenome->new(
	'xmlFiles'=>[
		'data/blastout1.xml',
		'data/blastout2.xml',
		'data/blastout3.xml',
		'data/blastout4.xml',
		'data/blastout5.xml',
		'data/blastout6.xml',
		'data/blastout7.xml',
		'data/blastout8.xml',
		'data/blastout9.xml',
		'data/blastout10.xm',
		'data/blastout11.xm',
		'data/blastout12.xm',
		'data/blastout13.xm',
		'data/blastout14.xm',
		'data/blastout15.xm',
		'data/blastout16.xm',
		'data/blastout17.xm',
		'data/blastout18.xm',
		'data/blastout19.xm',
		'data/blastout20.xm'
	],
	'numberOfCores'=>4,
	'resultArraySize'=>5,
	'percentIdentityCutoff'=>85,
	'coreGenomeThreshold'=>2,
	'outputDirectory'=>'temp/',
	'muscleExecutable'=>'/usr/bin/muscle',
	'accessoryType'=>'binary',
	'queryFile'=>'data/singleQueryFile.fasta'
);

ok(defined $pg, 'Object exists');
ok($pg->isa('Modules::PanGenome::PanGenome'), 'and it is the proper Modules::PanGenome::PanGenome class');

$pg->_generateOrderedNamesArray();

my @expectedOrderArray=('dbj|AP010960','gb|AFQI0100','gb|AKLB0100','gb|AKLN0100');
is($pg->_orderedNames->[0],$expectedOrderArray[0], "orderedNames[0] is the expected dbj|AP010960");
is($pg->_orderedNames->[1],$expectedOrderArray[1], "orderedNames[1] is the expected gb|AFQI0100");
is($pg->_orderedNames->[2],$expectedOrderArray[2], "orderedNames[2] is the expected gb|AKLB0100");
is($pg->_orderedNames->[3],$expectedOrderArray[3], "orderedNames[3] is the expected gb|AKLN0100");

#create filehandle
my $xmlFH = IO::File->new( '<' . $pg->xmlFiles->[0]) or die "$!";
my $xmler = Modules::Alignment::BlastResultFactory->new($xmlFH);

my $result;
for(1..5){
	$result= $xmler->nextResult;
}	
my ($accRef,$coreRef,$resultNumber)=$pg->_processQueue([$result], 1);
my @accData = split('\t',$accRef->[0]);

is($accData[1],1,"Correct value of accessory data for dbj|AP010960");
is($accData[2],1,"Correct value of accessory data for gb|AFQI0100");
is($accData[3],1,"Correct value of accessory data for gb|AKLB0100");
is($accData[4],1,"Correct value of accessory data for gb|AKLN0100");

foreach my $arrayRef(@{$coreRef}){
	foreach my $line(@{$arrayRef}){
		print "$line\n";
	}
}



1;



