#!/usr/bin/env perl


# Command being timed: "Panseq/t/output.t"
# User time (seconds): 696.37
# System time (seconds): 32.72
# Percent of CPU this job got: 98%
# Elapsed (wall clock) time (h:mm:ss or m:ss): 12:17.05
# Average shared text size (kbytes): 0
# Average unshared data size (kbytes): 0
# Average stack size (kbytes): 0
# Average total size (kbytes): 0
# Maximum resident set size (kbytes): 236900
# Average resident set size (kbytes): 0
# Major (requiring I/O) page faults: 75
# Minor (reclaiming a frame) page faults: 13096120
# Voluntary context switches: 7528
# Involuntary context switches: 7030
# Swaps: 0
# File system inputs: 7696
# File system outputs: 835440
# Socket messages sent: 0
# Socket messages received: 0
# Signals delivered: 0
# Page size (bytes): 4096
# Exit status: 0


use strict;
use warnings;
use Test::More tests=>19;
use Test::Pretty;
use File::Path qw/remove_tree/;
use Digest::MD5;
use IO::File;
use File::Copy;
use File::Basename;
use Getopt::Long;

my $type = $ARGV[0] // 'genomes';
my $removeRun = $ARGV[1] // 1;

#program locations
my $blastDirectory = '/usr/bin/';
my $mummerDirectory = '/usr/bin/';
my $muscleExecutable = '/usr/bin/muscle';

 GetOptions ('blastDirectory:s' => \$blastDirectory,
             'mummerDirectory:s' => \$mummerDirectory,
             'muscleExecutable:s' => \$muscleExecutable);

#get script location via File::Basename
my $SCRIPT_LOCATION = dirname(__FILE__);
print "SCRIPT_LOCATION: $SCRIPT_LOCATION\n";

my %plasmidsConfig=(
	queryDirectory=>"$SCRIPT_LOCATION/data/plasmids/",
	baseDirectory=>"$SCRIPT_LOCATION/plasmids/",
	numberOfCores=>1,
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>500,
	novelRegionFinderMode=>"no_duplicates",
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>'500',
	percentIdentityCutoff=>'90',
	coreGenomeThreshold=>'2',
	runMode=>'pan',
    nameOrId=>'name',
	overwrite=>1
);

my %queryConfig=(
	queryDirectory=>"$SCRIPT_LOCATION/data/genomes/",
	queryFile=>"$SCRIPT_LOCATION/data/testfragments.fasta",
	baseDirectory=>"$SCRIPT_LOCATION/query/",
	numberOfCores=>1,
	nameOrId=>'name',
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>1,
	novelRegionFinderMode=>"no_duplicates",
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>0,
	percentIdentityCutoff=>90,
	coreGenomeThreshold=>2,
	runMode=>'pan',
	storeAlleles=>1,
	overwrite=>1
);

my %genomesConfig=(
	queryDirectory=>"$SCRIPT_LOCATION/data/genomes/",
	baseDirectory=>"$SCRIPT_LOCATION/genomes/",
	numberOfCores=>1,
    nameOrId=>'name',
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>"1000",
	novelRegionFinderMode=>"no_duplicates",
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>'1000',
	percentIdentityCutoff=>'90',
	coreGenomeThreshold=>'3',
	runMode=>'pan',
	overwrite=>1
);

my %md5Sum=(
	plasmidsCoreSnps=>'afc0ccd1d1f8212b05ae78edb986fadf',
	plasmidsPanGenome=>'55475bf1f04dde9a52a9235190af53ef',
	plasmidsBinaryTable=>'04f2367db87c989b26c27d61c6908fa6',
	plasmidsSnpTable=>'4a8126a105bab8e7eeb0bde17b39e8b6',
	plasmidsBinaryPhylip=>'61c0c9a22041385779a5742944075729',
	plasmidsSnpPhylip=>'192df76ef39eb53a46405d7f17c60dcf',
	genomesCoreSnps=>'1ba7505faa64d1528cc5dba2dc420d34',
	genomesPanGenome=>'91dab20c6e944f580b9bd0ca59c2432d',
	genomesBinaryTable=>'b12ce2866684554a73b394d2ab71135f',
	genomesSnpTable=>'d612977784292d1286260e5e15a7892d',
	genomesBinaryPhylip=>'31aac5a26e78c41f06513f7181d86e2f',
	genomesSnpPhylip=>'2fe33bbd9d66d857ec7c42498e94f134',
	queryCoreSnps=>'8cff24ffeef699c9c07e98c9f486b46e',
	queryPanGenome=>'b770d4b30d6417ed43050f1254184bcb',
	queryBinaryTable=>'6cbe3340030974ce52f228031a4b30de',
	querySnpTable=>'3b10e7b6925e12a26eefc0d05763f1a6',
	queryBinaryPhylip=>'09e40c2e1be9d06d76c478418459982e',
	querySnpPhylip=>'8aaa390906f3ffbafc09413c441c48be',
	queryAlleles=>'41bfcb9626577dd46a2e0d1f2459c55d'
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
	
	#remove the ID column for testing, as it changes every run
	_removeIDColumn("$SCRIPT_LOCATION/$test/");
	
	my $md5 = _getMD5("$SCRIPT_LOCATION/$test/");
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

sub _getFilesFromDirectory{
	my $directory = shift;
	
	opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;	    
    return @dir;
}

#need to have test output be consistent
sub _removeIDColumn{
	my $directory = shift;
	
    my @dir = _getFilesFromDirectory($directory);
    
    foreach my $file(@dir){
    	unless(
    		$file eq 'pan_genome.txt'
    		|| $file eq 'core_snps.txt'
            || $file eq 'snp_table.txt'
    	){
    		next;
    	}
    	
    	my $originalFileName = $directory . $file;
    	my $modFileName = $originalFileName . 'mod';
    	
    	my $tempFH = IO::File->new('<' . $originalFileName) or die "Could not open $originalFileName";
    	my $tempOut = IO::File->new('>'. $directory . $file . 'mod') or die "Could not create modded file $modFileName";
    	
    	while(my $line = $tempFH->getline){
    		my @la = split("\t",$line);
    		shift @la;
    		$tempOut->print(join("\t",@la));
    	}
    	$tempOut->close();
    	$tempFH->close();
    	
    	#with File::Copy
    	move($modFileName,$originalFileName);
    }
	
}


sub _getMD5{
	my $directory=shift;	
	
    my @dir = _getFilesFromDirectory($directory);
	
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
	remove_tree("$SCRIPT_LOCATION/$t");
	unlink "$SCRIPT_LOCATION/$t.batch";
}

sub _runPanseq{
	my $t=shift;
	
	my $systemLine="perl $SCRIPT_LOCATION/../lib/panseq.pl $SCRIPT_LOCATION/$t.batch";
	print "Systemline: $systemLine\n";
	system($systemLine);
}


sub _createBatchFile{
	my $batchFile=shift;
	my $name=shift;
	
	my $batchFH=IO::File->new('>' . "$SCRIPT_LOCATION/$name.batch") or die "Could not create test batch file $name.batch";
	foreach my $key(keys %{$batchFile}){
		$batchFH->print("$key\t$batchFile->{$key}\n");
	}
	$batchFH->close();
}
