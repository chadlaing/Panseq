#!/usr/bin/env perl

use warnings;
use strict;
use Module::Build;
use IO::File;

my $build=Module::Build->new(
	module_name => 'Panseq',
	dist_version => '2016',
	license => 'perl',
	dist_abstract => 'Panseq: pan-genomic sequence analysis',
	dist_author => 'Chad Laing <chadlaing@inoutbox.com>',
	requires =>{
		'Carp' => 0,
		'perl' => '5.10.0',
		'Parallel::ForkManager' => 0,
		'Log::Log4perl' => 0,
		'FindBin' => 0,
		'IO::File' => 0,
		'File::Copy' => 0,
		'Role::Tiny' => 0,
		'Bio::SeqIO' => 0,
		'Bio::DB::Fasta' => 0,
		'File::Temp' => 0,
		'Bio::Seq' => 0,
		'File::Basename'=>0,
		'Tie::Log4perl'=>0,
		'Archive::Zip'=>0,
		'Test2::Bundle::Extended'=>0,
		'Getopt::Long'=>0,
		'Digest::MD5'=>0
	}
);
$build->create_build_script();

