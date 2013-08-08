#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Log::Log4perl;
use Modules::NovelRegion::NovelRegionFinder;
Log::Log4perl->init('simple.conf');

my $logger =  Log::Log4perl->get_logger();

my $nrf = Modules::NovelRegion::NovelRegionFinder->new(
		mode=>'unique',
		coordsFile=>'/home/chad/panseq/output/ecoli_test/unique_regions.coords',
		queryFile=>'/home/chad/panseq/output/ecoli_test/lastNovelRegionsFile.fasta',
		minimumNovelRegionSize=>1000
	);
	$nrf->findNovelRegions();