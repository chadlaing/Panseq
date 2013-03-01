#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Modules::Setup::Panseq;
use Modules::Setup::Settings;
use Log::Log4perl qw/:easy/;

#need to get the settings up-front to allow log-file output
my $settings = Modules::Setup::Settings->new($ARGV[0]);

#MUMmmer defaults its messages to STDERR
#we want them logged
#closing STDERR and associating it with Log4perl is done below
# close STDERR;
# tie *STDERR, "Tie::Log4perl";

my $panseq = Modules::Setup::Panseq->new($settings);
Log::Log4perl->init("$FindBin::Bin/log4p.conf");

$panseq->run();


#HOOKS for log4p.conf
sub nucmer_run_file{
	return ($settings->baseDirectory . 'logs/NucmerRun.pm.log'); 
}

sub panseq_files_file{
	return ($settings->baseDirectory . 'logs/PanseqFiles.pm.log'); 
}

sub novel_region_finder_file{
	return ($settings->baseDirectory . 'logs/NovelRegionFinder.pm.log'); 
}

sub sequence_name_file{
	return ($settings->baseDirectory . 'logs/SequenceName.pm.log'); 
}

sub multi_fasta_sequence_name_file{
	return ($settings->baseDirectory . 'logs/MultiSequenceFastaName.pm.log'); 
}

sub novel_iterator_file{
	return ($settings->baseDirectory . 'logs/NovelIterator.pm.log'); 
}



