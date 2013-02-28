#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Test::More tests=>15;
use Test::Pretty;

use Modules::Fasta::SequenceName;

#test the SequenceName.pm module
my $sn = Modules::Fasta::SequenceName->new('gi|148566106|gb|CP000711.1| Enterobacteria phage CUS-3, complete genome');

ok(defined $sn, 'Object exists');
ok($sn->isa('Modules::Fasta::SequenceName'), 'and it is the proper Modules::Fasta::SequenceName class');

#test it requires a single argument
my $error = 'constructor requires a single value';
eval {my $obj = Modules::Fasta::SequenceName->new()};
is($@ =~ m/$error/, 1, 'Properly caught no argument in Modules::Fasta::SequenceName constructor');

eval {my $obj = Modules::Fasta::SequenceName->new('single argument')};
is($@ =~ m/$error/,'', 'Properly initialized with one argument in Modules::Fasta::SequenceName');

eval {my $obj = Modules::Fasta::SequenceName->new('first','second')};
is($@ =~ m/$error/,1, 'Properly caught multiple argument in Modules::Fasta::SequenceName constructor');

#test the various ways the data should be handled
#gi|148566106|gb|CP000711.1| Enterobacteria phage CUS-3, complete genome
my %testValues=(
	'name'=>{
			'string'=>'name=|test1|gi|148566106|gb|CP000711.1| Enterobacteria phage CUS-3, complete genome',
			'correct'=>'test1'
	},
	'lcl'=>{
			'string'=>'lcl|test1|gi|148566106|gb|CP000711.1| Enterobacteria phage CUS-3, complete genome',
			'correct'=>'test1'
	},
	'ref'=>{
			'string'=>'ref|NC_123456|gi|148566106|gb|CP000711.1| Enterobacteria phage CUS-3, complete genome',
			'correct'=>'ref|NC_123456'
	},
	'gb'=>{
			'string'=>'gi|148566106|gb|CP000711.1| Enterobacteria phage CUS-3, complete genome',
			'correct'=>'gb|CP000711'
	},
	'emb'=>{
			'string'=>'gi|148566106|emb|CP000711.1| Enterobacteria phage CUS-3, complete genome',
			'correct'=>'emb|CP000711'
	},
	'dbj'=>{
			'string'=>'gi|148566106|dbj|CP000711.1| Enterobacteria phage CUS-3, complete genome',
			'correct'=>'dbj|CP000711'
	},
	'Segment'=>{
			'string'=>'gi_148566106_dbj_CP000711.1| Enterobacteria phage CUS-3, complete genome|Segment=125',
			'correct'=>'gi_148566106_dbj_CP000711.1| Enterobacteria phage CUS-3, complete genome'
	},
	'Length'=>{
			'string'=>'gi_148566106_dbj_CP000711.1| Enterobacteria phage CUS-3, complete genome|Length=1245555',
			'correct'=>'gi_148566106_dbj_CP000711.1| Enterobacteria phage CUS-3, complete genome'
	},
	'original'=>{
			'string'=>'gi_148566106_dbj_CP000711_1_ Enterobacteria phage CUS-3, complete genome_Length=1245555',
			'correct'=>'gi_148566106_dbj_CP000711_1_ Enterobacteria phage CUS-3, complete genome_Length=1245555'
	},
	'gi'=>{			
		'string'=>'gi|148566106|gret|CP000711.1| Enterobacteria phage CUS-3, complete genome',
		'correct'=>'gi|148566106'
	}
);


foreach my $test(keys %testValues){
	my $sn = Modules::Fasta::SequenceName->new($testValues{$test}->{'string'});
	is($sn->name,$testValues{$test}->{'correct'}, 'Name is the expected ' . $testValues{$test}->{'correct'} . ' in Modules::Fasta::SequenceName->name(' . $test . ')');
}

