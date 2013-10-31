#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Test::More tests=>19;
use Test::Pretty;
use File::Path qw/remove_tree/;
use Digest::MD5;
use IO::File;

my $type = $ARGV[0] // 'genomes';
my $removeRun = $ARGV[1] // 1;

#program locations
my $blastDirectory = '/usr/bin/';
my $mummerDirectory = '/usr/bin/';
my $muscleExecutable = '/usr/bin/muscle';

my %plasmidsConfig=(
	queryDirectory=>"$FindBin::Bin/data/plasmids/",
	baseDirectory=>"$FindBin::Bin/plasmids/",
	numberOfCores=>4,
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>"500",
	novelRegionFinderMode=>"no_duplicates",
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>'500',
	percentIdentityCutoff=>'90',
	coreGenomeThreshold=>'3',
	runMode=>'pan'
);

my %queryConfig=(
	queryDirectory=>"$FindBin::Bin/data/genomes/",
	queryFile=>"$FindBin::Bin/data/testfragments.fasta",
	baseDirectory=>"$FindBin::Bin/query/",
	numberOfCores=>"1",
	nameOrId=>'name',
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>"500",
	novelRegionFinderMode=>"no_duplicates",
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>'0',
	percentIdentityCutoff=>'90',
	coreGenomeThreshold=>'2',
	runMode=>'pan',
	storeAlleles=>1
);

my %genomesConfig=(
	queryDirectory=>"$FindBin::Bin/data/genomes/",
	baseDirectory=>"$FindBin::Bin/genomes/",
	numberOfCores=>"1",
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>"1000",
	novelRegionFinderMode=>"no_duplicates",
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>'1000',
	percentIdentityCutoff=>'90',
	coreGenomeThreshold=>'3',
	runMode=>'pan'
);

my %md5Sum=(
	plasmidsCoreSnps=>'8c167d9429ec18aca2c84190d97e346a',
	plasmidsPanGenome=>'836941e2bb89f6540a36d04610fff1c9',
	plasmidsBinaryTable=>'a42ac19fd24b06d92cda1a9cfea8112a',
	plasmidsSnpTable=>'e29f54adc12c207b523f8a8c02109b1a',
	plasmidsBinaryPhylip=>'a033883cc4bcf457fefca7da76d942e6',
	plasmidsSnpPhylip=>'eedf538b1f01a8afa53bc813c26eb03b',
	genomesCoreSnps=>'4c25c6ba2475f5d9cc6641131dfd05c3',
	genomesPanGenome=>'478bbf41ec3a894b5a0a289ee6d4c132',
	genomesBinaryTable=>'c1002e14a5da5dad544ce53d8e93a2e2',
	genomesSnpTable=>'2ec5fe0256aa273c46c69a6ca7364c76',
	genomesBinaryPhylip=>'2e126c05a39fc00eb54c5a74fb2879fc',
	genomesSnpPhylip=>'e0041d87a497056b7934b9fa181820a7',
	queryCoreSnps=>'665e547de347114db35fa8a813210956',
	queryPanGenome=>'7d429141afd723524380342c293f6a62',
	queryBinaryTable=>'7b7505fee1407c86f5b8f75d03c9b24c',
	querySnpTable=>'d41d8cd98f00b204e9800998ecf8427e',
	queryBinaryPhylip=>'14b852f4f00745ae8e8190a8fdb9dbe2',
	querySnpPhylip=>'d41d8cd98f00b204e9800998ecf8427e',
	queryAlleles=>'201833090475055cd4ec9f28ce4adc70'
);

#create the Batch files and test the output of Panseq to ensure no breaking changes have occurred
#generate data first, so all tests are at the bottom of the output
foreach my $test(@{['plasmids','query','genomes']}){
	if($type ne 'genomes' and $test eq 'genomes'){
		next;
	}
	my %config;
	if($test eq 'genomes'){
		%config = %genomesConfig;
	}
	elsif($test eq 'query'){
		%config = %queryConfig;
	}
	elsif($test eq 'plasmids'){
		%config = %plasmidsConfig;
	}	
	
	_createBatchFile(\%config,$test);
	_runPanseq($test);
}

#compare the digests of the files for correctness
foreach my $test(@{['plasmids','query','genomes']}){
	if($type ne 'genomes' and $test eq 'genomes'){
		next;
	}
	
	my $md5 = _getMD5("$FindBin::Bin/$test/");
	is($md5->{'coreSnps'},$md5Sum{"${test}CoreSnps"},"${test}CoreSnps generated correctly");
	is($md5->{'panGenome'},$md5Sum{"${test}PanGenome"},"${test}PanGenome generated correctly");
	is($md5->{'binaryTable'},$md5Sum{"${test}BinaryTable"},"${test}BinaryTable generated correctly");
	is($md5->{'snpTable'},$md5Sum{"${test}SnpTable"},"${test}SnpTable generated correctly");
	is($md5->{'snpPhylip'},$md5Sum{"${test}SnpPhylip"},"${test}SnpPhylip generated correctly");
	is($md5->{'binaryPhylip'},$md5Sum{"${test}BinaryPhylip"},"${test}BinaryPhylip generated correctly");
	
	if($test eq 'query'){
		is($md5->{'locusAlleles'},$md5Sum{"${test}Alleles"},"${test}Alleles generated correctly");
	}
	
	if($removeRun == 1){
		_removeRun($test);
	}
	
}

sub _getMD5{
	my $directory=shift;
	
	opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;	
	
	my %md5Hash;
	my $digester = Digest::MD5->new();
    foreach my $fileName(sort @dir){
    	
        if((substr( $fileName, 0, 1 ) eq '.') || (-d $directory . $fileName)){
        	next;
        }
        
    	my $inFH=IO::File->new('<' . $directory . $fileName) or die "Could not open $directory$fileName";
    	$inFH->binmode();   	
 
        my $md5sum = $digester->addfile($inFH)->hexdigest;
     
        if($fileName eq 'core_snps.txt'){
        	$md5Hash{'coreSnps'}=$md5sum;
        }
        elsif($fileName eq 'pan_genome.txt'){
        	$md5Hash{'panGenome'}=$md5sum;
        }
        elsif($fileName eq 'binary_table.txt'){
        	$md5Hash{'binaryTable'}=$md5sum;
        }
        elsif($fileName eq 'snp_table.txt'){
        	$md5Hash{'snpTable'}=$md5sum;
        }
        elsif($fileName eq 'binary.phylip'){
        	$md5Hash{'binaryPhylip'}=$md5sum;
        }
        elsif($fileName eq 'snp.phylip'){
        	$md5Hash{'snpPhylip'}=$md5sum;
        }
        elsif($fileName eq 'locus_alleles.fasta'){
        	$md5Hash{'locusAlleles'}=$md5sum;
        }        
        $inFH->close();
    }
    return \%md5Hash;	
}

sub _removeRun{
	my $t=shift;
	
	#with File::Path
	remove_tree("$FindBin::Bin/$t");
	unlink "$FindBin::Bin/$t.batch";
}

sub _runPanseq{
	my $t=shift;
	
	my $systemLine="perl $FindBin::Bin/../lib/panseq.pl $t.batch";
	system($systemLine);
}


sub _createBatchFile{
	my $batchFile=shift;
	my $name=shift;
	
	my $batchFH=IO::File->new('>' . "$FindBin::Bin/$name.batch") or die "Could not create test batch file $name.batch";
	foreach my $key(keys %{$batchFile}){
		$batchFH->print("$key\t$batchFile->{$key}\n");
	}
	$batchFH->close();
}
