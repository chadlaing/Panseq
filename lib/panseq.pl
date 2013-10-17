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
close STDERR;
tie *STDERR, "Tie::Log4perl";

my $panseq = Modules::Setup::Panseq->new($settings);
Log::Log4perl->init("$FindBin::Bin/log4p.conf");

$panseq->run();


#HOOKS for log4p.conf
sub panseq_master_log_file{
	return ($settings->baseDirectory . 'Master.log');
}



