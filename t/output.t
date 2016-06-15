#!/usr/bin/env perl

use strict;
use warnings;
#use Test::More tests=>26;
#use Test::Pretty;
use Test2::Bundle::Extended;
use File::Path qw/remove_tree/;
use Digest::MD5;
use IO::File;
use File::Copy;
use File::Basename;
use Getopt::Long;

#set up plan
plan 26;

#test options
my $blastDirectory = '/usr/bin/';
my $mummerDirectory = '/usr/bin/';
my $muscleExecutable = '/usr/bin/muscle';
my $numberOfCores = 1;
my $type = 'genomes';
my $removeRun = 1;

 GetOptions ('blastDirectory=s' => \$blastDirectory,
             'mummerDirectory=s' => \$mummerDirectory,
             'muscleExecutable=s' => \$muscleExecutable,
             'type=s' => \$type,
             'removeRun=i' => \$removeRun,
             #Because of the change in splitting the fasta file, the number of cores
             #affects the order of fasta sequences, which is guaranteed to be the same for each run
             #using the same number of cores, but not the same among different core numbers
             #therefore, for consistent testing, a single core is used
             #'numberOfCores=i' => \$numberOfCores
             );

#get script location via File::Basename
my $SCRIPT_LOCATION = dirname(__FILE__);
print "SCRIPT_LOCATION: $SCRIPT_LOCATION\n";

my %plasmidsConfig=(
	queryDirectory=>"$SCRIPT_LOCATION/data/plasmids/",
	baseDirectory=>"$SCRIPT_LOCATION/plasmids/",
	numberOfCores=>$numberOfCores,
	mummerDirectory=>$mummerDirectory,
	blastDirectory=>$blastDirectory,
	minimumNovelRegionSize=>500,
	muscleExecutable=>$muscleExecutable,
	fragmentationSize=>'500',
	percentIdentityCutoff=>'90',
	coreGenomeThreshold=>'2',
	runMode=>'pan',
    nameOrId=>'name',
	overwrite=>1,
    storeAlleles=>1,
    allelesToKeep=>2
);

my %queryConfig=(
	queryDirectory=>"$SCRIPT_LOCATION/data/genomes/",
	queryFile=>"$SCRIPT_LOCATION/data/testfragments.fasta",
	baseDirectory=>"$SCRIPT_LOCATION/query/",
	numberOfCores=>$numberOfCores,
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
	numberOfCores=>$numberOfCores,
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
	plasmidsCoreSnps=>'a7d2902d80543446a6701e8f8770301b',
    plasmidsPanGenome=>'168d75b59dbe825cd91222906a4f5645',
    plasmidsBinaryTable=>'fbfffe4e58a1dfc1cd04cb29b5146c0d',
    plasmidsSnpTable=>'188356ccf7f73066c333e82dd18e531d',
    plasmidsBinaryPhylip=>'db5e15a38d9b7b7be53811df302d7558',
    plasmidsSnpPhylip=>'6c8f15448b0f19be9efbb79e16256350',
    plasmidsCoreFragments=>'2fa21523c2c0e9abde0836f2a754640e',
    plasmidsAccessoryFragments=>'f17e29fd8ca3dbaac3033ce188018465',
    plasmidsNameConversion=>'da9678fa95a0def763ad014ec7153779',
	genomesCoreSnps=>'83beb7f1fbdcf2f6cb38cb5604b8385b',
    genomesPanGenome=>'3e00bb9e7d7fa9b02b34052fd005fa00',
    genomesBinaryTable=>'1f1aaef9c674a5e847cae718964b0385',
    genomesSnpTable=>'63ca2aa391938fb521375f0b2353bb06',
    genomesBinaryPhylip=>'4b341c515a3aa54377b7a7f8a9e71d17',
    genomesSnpPhylip=>'f438ce263872f10b36cb3247a3d5dd59',
    genomesCoreFragments=>'117d52a380e05eddd33a31d07a4f7829',
    genomesAccessoryFragments=>'6dca4cb62aabfbca4d54279d959fc451',
    genomesNameConversion=>'e90cc17adc92f2d63106d58dff86860a',
	queryCoreSnps=>'67daf16fc332d162e485227f80e22958',
    queryPanGenome=>'a603a5526709da34bc854363045c94bf',
    queryBinaryTable=>'1727cd2ef07eb6082793717521d7146f',
    querySnpTable=>'aca086386f540630499b24f0e3644e3d',
    queryBinaryPhylip=>'183fea98a21e4f9eae54e486f1f08821',
    querySnpPhylip=>'972f68dff31780ddab9c1ac87b78eddd',
    queryAlleles=>'7aec36d7ee53447e0dd5e82be3d2f9bc',
    queryNameConversion=>'e90cc17adc92f2d63106d58dff86860a'
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

    #remove fast headers for Fragments files, as they include the IDs, which change every run

	my $md5 = _getMD5("$SCRIPT_LOCATION/$test/");
	is($md5->{'coreSnps'},$md5Sum{"${test}CoreSnps"},"${test}CoreSnps generated correctly");
	is($md5->{'panGenome'},$md5Sum{"${test}PanGenome"},"${test}PanGenome generated correctly");
	is($md5->{'binaryTable'},$md5Sum{"${test}BinaryTable"},"${test}BinaryTable generated correctly");
	is($md5->{'snpTable'},$md5Sum{"${test}SnpTable"},"${test}SnpTable generated correctly");
	is($md5->{'snpPhylip'},$md5Sum{"${test}SnpPhylip"},"${test}SnpPhylip generated correctly");
	is($md5->{'binaryPhylip'},$md5Sum{"${test}BinaryPhylip"},"${test}BinaryPhylip generated correctly");
    is($md5->{'nameConversion'},$md5Sum{"${test}NameConversion"},"${test}NameConversion generated correctly");
	
	if($test eq 'query'){
		is($md5->{'locusAlleles'},$md5Sum{"${test}Alleles"},"${test}Alleles generated correctly");
	}
    else{
        is($md5->{'accessoryFragments'},$md5Sum{"${test}AccessoryFragments"},"${test}AccessoryFragments generated correctly");
        is($md5->{'coreFragments'},$md5Sum{"${test}CoreFragments"},"${test}CoreFragments generated correctly");
    }


	if($removeRun == 1){
		_removeRun($test);
	}
	
}

#end
done_testing;

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
            || $file eq 'accessoryGenomeFragments.fasta'
            || $file eq 'coreGenomeFragments.fasta'
    	){
    		next;
    	}
    	
    	my $originalFileName = $directory . $file;
    	my $modFileName = $originalFileName . 'mod';
    	
    	my $tempFH = IO::File->new('<' . $originalFileName) or die "Could not open $originalFileName";
    	my $tempOut = IO::File->new('>'. $directory . $file . 'mod') or die "Could not create modded file $modFileName";
    	
        if($file eq 'accessoryGenomeFragments.fasta' || $file eq 'coreGenomeFragments.fasta'){
            while(my $line = $tempFH->getline){
                if($line =~ m/^>/){
                    next;
                }
                else{
                    $tempOut->print($line);
                }
            }
        }
        else{
            while(my $line = $tempFH->getline){
                my @la = split("\t",$line);
                shift @la;
                $tempOut->print(join("\t",@la));
            }
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
        elsif($fileName eq 'accessoryGenomeFragments.fasta'){
            $md5Hash{'accessoryFragments'}=$md5sum;
        }
        elsif($fileName eq 'coreGenomeFragments.fasta'){
            $md5Hash{'coreFragments'}=$md5sum;
        }
        elsif($fileName eq 'phylip_name_conversion.txt'){
            $md5Hash{'nameConversion'}=$md5sum;
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
