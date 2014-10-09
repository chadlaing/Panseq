#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Test::More tests=>7;
use Test::Pretty;
use Log::Log4perl qw/:easy/;
use Modules::NovelRegion::NovelRegionFinder;
use Modules::Setup::Settings;

my $settings = Modules::Setup::Settings->new();


Log::Log4perl->easy_init($DEBUG);
my $logger =  Log::Log4perl->get_logger();
$logger->info("Testing");
#data


#test the NovelRegionFinder.pm module
my $nrf = Modules::NovelRegion::NovelRegionFinder->new();
#add the required setting elements, as we are not actually loading a file
$nrf->settings($settings);
$nrf->settings->_genomeNameFromContig({
    query=>'query',
    query2=>'query2'
    });



ok(defined $nrf, 'Object exists');
ok($nrf->isa('Modules::NovelRegion::NovelRegionFinder'), 'and it is the proper Modules::NovelRegion::NovelRegionFinder class');

#_updateComparsionHash
#query:   1.................................................................5000
#ref:           151............2002               3006...............4817
$nrf->_updateComparisonHash('151	2002	151	2002	1852	1852	100	6000	5000	ref	query');
is($nrf->_comparisonHash->{'query'}->{'ref'},',151..2002','_updateComparsionHash correctly added a new sequence');

$nrf->_updateComparisonHash('3006	4817	3006	4817	1810	1810	100	6000	5000	ref	query');
is($nrf->_comparisonHash->{'query'}->{'ref'},',151..2002,3006..4817','_updateComparsionHash correctly added to an existing sequence');

#create a _queryFastaHeadersHash for testing
$nrf->_queryFastaHeadersHash->{'query'}=5000;
#check for an empty hit
$nrf->_queryFastaHeadersHash->{'query2'}=4500;

my $novelRegions = $nrf->_getNoDuplicates();
is($novelRegions->{'query'},',1..150,2003..3005,4818..5000', ',1..150,2003..3005,4818..5000 correctly identified as novelRegions');
is($novelRegions->{'query2'},',1..4500','Query2 with no hit, correctly identified as 1..4500');


#test for query vs multiple ref name hits
#ditto the settings setup
my $nrf2 = Modules::NovelRegion::NovelRegionFinder->new();
$nrf2->settings($settings);
$nrf2->settings->_genomeNameFromContig({
    query=>'query',
    query2=>'query2'
    });
#query:     1.............................................................................7000
#ref1:      1..........1000
#ref2:              500................2001
#ref3:							1500........................4320    5000....5500
#ref4:															                  5630..6900
$nrf2->_queryFastaHeadersHash->{'query'}=7000;

$nrf2->_updateComparisonHash('1	1000	1	1000	1000	1000	100	1000	7000	ref1	query');
$nrf2->_updateComparisonHash('1	1501	500	2001	1501	1501	100	1501	7000	ref2	query');
$nrf2->_updateComparisonHash('1	2821	1500	4320	2821	2821	100	6000	7000	ref3	query');
$nrf2->_updateComparisonHash('3000	3500	5000	5500	501	501	100	6000	7000	ref3	query');
$nrf2->_updateComparisonHash('100	1370	5630	6900	1271	1271	100	2000	7000	ref4	query');
my $novelRegions2 = $nrf2->_getNoDuplicates();
is($novelRegions2->{'query'},',4321..4999,5501..5629,6901..7000', ',4321..4999,5501..5629,6901..7000 correctly identified as novelRegions from multiple refs');

$logger->info("Finished testing NovelRegionFinder");





