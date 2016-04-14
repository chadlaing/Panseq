#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Net::FTP;
use File::Basename;

#get script location via File::Basename
my $SCRIPT_LOCATION = dirname(__FILE__);

my $localFile = $ARGV[0] // $SCRIPT_LOCATION . '/assembly_summary.txt';

#sets up the parameters for ncbi ftp connection
my $host = 'ftp.ncbi.nlm.nih.gov';

#new FTP directory structure
my $genomeSummaryFile = '/genomes/refseq/bacteria/assembly_summary.txt';

#constructs the connection
my $ftp = Net::FTP->new($host, Debug => 1,Passive => 1, Timeout => 1000) or die "Cannot connect to genbank: $@";
#log in as anonymous, use email as password
$ftp->login("anonymous",'chadlaing@inoutbox.com') or die "Cannot login ", $ftp->message;

$ftp->binary();
$ftp->get($genomeSummaryFile, $localFile) or die ("download of $genomeSummaryFile failed\n");
$ftp->ascii();
