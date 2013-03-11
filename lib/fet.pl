#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin";
use warnings;
use strict;
use Log::Log4perl;
Log::Log4perl->easy_init();


use Modules::Stats::FET;

my $fet = Modules::Stats::FET->new(
	'testCharacters'=>['A','T','C','G'],
	'inputFile'=>'/home/chad/workspace/segenomes/core_snps_edit.txt',
	'excludedCharacters'=>['-'],
	'group1'=>['SE_9-698','SE_11-355'
	],
	'group2'=>[
		'SE_10-101','SE_11-1095','SE_11-1174','SE_11-1175','SE_11-354','SE_11-356','SE_11-360','SE_11-361',
		'SE_11-221','SE_11-353','SE_11-357','SE_11-358','SE_11-359','SE_12-8','SE_12-9','SE_9-641'
	]
);
$fet->run();