#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Modules::Setup::Panseq;
use Log::Log4perl qw/:easy/;
Log::Log4perl->easy_init($INFO);

my $panseq = Modules::Setup::Panseq->new($ARGV[0]);
$panseq->run();


# #HOOKS for log4p.conf
# sub nucmer_run_file{
# 	return ($settings->baseDirectory . 'logs/NucmerRun.pm.log'); 
# }

# sub panseq_files_file{
# 	return ($settings->baseDirectory . 'logs/PanseqFiles.pm.log'); 
# }

# sub novel_region_finder_file{
# 	return ($settings->baseDirectory . 'logs/NovelRegionFinder.pm.log'); 
# }

# sub sequence_name_file{
# 	return ($settings->baseDirectory . 'logs/SequenceName.pm.log'); 
# }

# sub multi_fasta_sequence_name_file{
# 	return ($settings->baseDirectory . 'logs/MultiSequenceFastaName.pm.log'); 
# }

# sub novel_iterator_file{
# 	return ($settings->baseDirectory . 'logs/NovelIterator.pm.log'); 
# }



