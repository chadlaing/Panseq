#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use Modules::Setup::Panseq;
use Modules::Setup::Settings;
use Log::Log4perl qw/:easy/;

#need to get the settings up-front to allow log-file output
my $settings = Modules::Setup::Settings->new($ARGV[0]);
print "settings file: $ARGV[0]\n";

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

sub blast_run_file{
	return ($settings->baseDirectory . 'logs/BlastRun.pm.log'); 
}

sub blast_result_factory_file{
	return ($settings->baseDirectory . 'logs/BlastResultFactory.pm.log'); 
}

sub snp_finder_file{
	return ($settings->baseDirectory . 'logs/SNPFinder.pm.log'); 
}

sub fasta_file_splitter_file{
	return ($settings->baseDirectory . 'logs/FastaFileSplitter.pm.log'); 
}

sub sequence_retriever_file{
	return ($settings->baseDirectory . 'logs/SequenceRetriever.pm.log'); 
}

sub segment_maker_file{
	return ($settings->baseDirectory . 'logs/SegmentMaker.pm.log'); 
}

sub pan_genome_file{
	return ($settings->baseDirectory . 'logs/PanGenome.pm.log'); 
}



