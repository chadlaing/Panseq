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
use File::Copy;

my $type = $ARGV[0] // 'genomes';
my $removeRun = $ARGV[1] // 1;

#program locations
my $blastDirectory = '/usr/bin/';
my $mummerDirectory = '/usr/bin/';
my $muscleExecutable = '/usr/bin/muscle';

my %plasmidsConfig=(
	queryDirectory=>"$FindBin::Bin/data/plasmids/",
	baseDirectory=>"$FindBin::Bin/plasmids/",
	numberOfCores=>1,
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
	numberOfCores=>1,
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
	numberOfCores=>1,
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
	plasmidsCoreSnps=>'225996c682a42630917f7c9e917bfc30',
	plasmidsPanGenome=>'6c398139935d400379728728eb15de95',
	plasmidsBinaryTable=>'b26574a7251f2fc0060e65cba1dbbc53',
	plasmidsSnpTable=>'3a88ed79ddcf0bc525319de29f4b8e23',
	plasmidsBinaryPhylip=>'a033883cc4bcf457fefca7da76d942e6',
	plasmidsSnpPhylip=>'eedf538b1f01a8afa53bc813c26eb03b',
	genomesCoreSnps=>'d0c99a6ea5d35bbd5d275d9a7fdd1b1b',
	genomesPanGenome=>'44619ab2f10d75509b8a8993cb900223',
	genomesBinaryTable=>'cfb2368af43146c273109c1f779d44db',
	genomesSnpTable=>'25343dcc1e33fd4bbc0195542ed35a6b',
	genomesBinaryPhylip=>'e576f7865893bf8e00bea0081b5f1aea',
	genomesSnpPhylip=>'14dd6176feba701c23589b68c025bd2c',
	queryCoreSnps=>'c0df4f4122d388a6d498f61910060436',
	queryPanGenome=>'be84039b4c11f7e5e731105d6fdf3378',
	queryBinaryTable=>'6d2e38b5ab20c77ad57f7ea02b7cf26a',
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
	
	#remove the ID column for testing, as it changes every run
	_removeIDColumn("$FindBin::Bin/$test/");
	
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

sub _getFilesFromDirectory{
	my $directory = shift;
	
	opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;	    
    return @dir;
}

#need to have test output ber consistent
sub _removeIDColumn{
	my $directory = shift;
	
    my @dir = _getFilesFromDirectory($directory);
    
    foreach my $file(@dir){
    	unless(
    		$file eq 'pan_genome.txt'
    		|| $file eq 'core_snps.txt'
    		|| $file eq 'binary_table.txt'
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
