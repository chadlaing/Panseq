#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/";
use Test::More tests=>7;
use Test::Pretty;
use File::Path qw/remove_tree/;
use Digest::MD5;
use IO::File;

my $type = $ARGV[0] // 'genomes';

#program locations
my $blastDirectory = '/usr/bin/';
my $mummerDirectory = '/usr/bin/';
my $muscleExecutable = '/usr/bin/muscle';

my %plasmidsConfig=(
	queryDirectory=>"$FindBin::Bin/data/plasmids/",
	baseDirectory=>"$FindBin::Bin/plasmids/",
	numberOfCores=>"4",
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
	numberOfCores=>"4",
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>"500",
	novelRegionFinderMode=>"no_duplicates",
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>'0',
	percentIdentityCutoff=>'90',
	coreGenomeThreshold=>'3',
	runMode=>'pan'
);

my %genomesConfig=(
	queryDirectory=>"$FindBin::Bin/data/genomes/",
	baseDirectory=>"$FindBin::Bin/genomes/",
	numberOfCores=>"4",
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
	plasmidsCoreSnps=>'dd7231bfcf77ea3c7883bbfb3818d4fe',
	plasmidsPanGenome=>'61742e5586810d22ce66fc7b9587a8a7',
	plasmidsBinaryTable=>'2f686ed44db1809e22aa2e76f6bce984',
	plasmidsSnpTable=>'ae80bc872dcde45839329e7bd165145e',
	plasmidsBinaryPhylip=>'ba4bb0049b6d592755cd52853d05f89c',
	plasmidsSnpPhylip=>'a9a66896db1800f9cfafafa2913804e6',
	genomesCoreSnps=>'30c1123981af005875348fee9223c402',
	genomesPanGenome=>'894765f4508cea9f8358da582f93ef65',
	genomesBinaryTable=>'30a955b59b67ee7a8550169ae56e0f51',
	genomesSnpTable=>'b2dcbf4a6e09b5ba870e678556228560',
	genomesBinaryPhylip=>'f7dd22ffc7736219298f668933b6caa6',
	genomesSnpPhylip=>'6c8524ec1478aa5a929ded90096ee99f',
	queryCoreSnps=>'fbf65d84b2416ee35629f2d53f67fc8a',
	queryPanGenome=>'1a2cbde5f3ef098c3ccf4daae35bad0d',
	queryBinaryTable=>'f7d06debab36e48c823de151d6a62401',
	querySnpTable=>'d41d8cd98f00b204e9800998ecf8427e',
	queryBinaryPhylip=>'75058ee0a94f9edbd20c7fe1cb2f0f8c',
	querySnpPhylip=>'d41d8cd98f00b204e9800998ecf8427e',
	queryAlleles=>'7cf30ab8f2368a37b26672f8de5ecb82'
);


#create the Batch files and test the output of Panseq to ensure no breaking changes have occurred
_createBatchFile(\%plasmidsConfig,'plasmids');
_runPanseq('plasmids');

_createBatchFile(\%queryConfig,'query');
_runPanseq('query');


if($type eq 'genomes'){
	_createBatchFile(\%genomesConfig,'genomes');
	_runPanseq('genomes');
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
}
_removeRun('query');
_removeRun('plasmids');

if($type eq 'genomes'){
	_removeRun('genomes');
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
      	print "File: $fileName md5sum: $md5sum\n";
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
