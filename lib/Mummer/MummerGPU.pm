#!/usr/bin/perl

#takes in the queryFile (all query files in one file)
#and the referenceFile (all reference files in one file)
#mummer needs the single query file with all seqs vs each ref sequence in a separate file

package Mummer::MummerGPU;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../";
use Carp;
use IO::File;
use File::Temp;
use Parallel::ForkManager;
use FileInteraction::FileManipulation;
use FileInteraction::Fasta::SequenceRetriever;
use Log::Log4perl;
use Bio::SeqIO;
use File::Copy;

#object creation
sub new{
	my $class=shift;
    my $self = {};
    bless ($self, $class);
    $self->_initialize(@_);   
    return $self;
}

#class variables

sub _queryFile{
	my $self=shift;
	$self->{'_queryFile'}=shift // return $self->{'_queryFile'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub _numberOfCores{
	my $self=shift;
	$self->{'_numberOfCores'}=shift // return $self->{'_numberOfCores'};
}

sub deltaFile{
	my $self=shift;
	$self->{'_deltaFile'}=shift // return $self->{'_deltaFile'};
}

sub _mummerDirectory{
	my $self=shift;
	$self->{'_mummerDirectory'}=shift // return $self->{'_mummerDirectory'};
}

sub _baseDirectory{
	my $self=shift;
	$self->{'_baseDirectory'}=shift // return $self->{'_baseDirectory'};
}

sub _referenceFileArray{
	my $self=shift;
	$self->{'_referenceFile'}=shift // return $self->{'_referenceFile'};
}

sub _refFileSizeLimit{
	my $self=shift;
	$self->{'_refFileSizeLimit'}=shift // return $self->{'_refFileSizeLimit'};	
}

sub _systemLineBase{
	my $self=shift;
	$self->{'_systemLineBase'}=shift // return $self->{'_systemLineBase'};
}

sub _gpu{
	my $self=shift;
	$self->{'_gpu'}=shift // return $self->{'_gpu'};
}


sub _l{
	my $self=shift;
	$self->{'_l'}=shift // return $self->{'_l'};
}

sub _g{
	my $self=shift;
	$self->{'_g'}=shift // return $self->{'_g'};
}

sub _d{
	my $self=shift;
	$self->{'_d'}=shift // return $self->{'_d'};
}

sub _c{
	my $self=shift;
	$self->{'_c'}=shift // return $self->{'_c'};
}

sub _b{
	my $self=shift;
	$self->{'_b'}=shift // return $self->{'_b'};
}

sub _p{
	my $self=shift;
	$self->{'_p'}=shift // return $self->{'_p'};
}

sub coordsFile{
	my $self=shift;
	$self->{'_coordsFile'}=shift // return $self->{'_coordsFile'};
}

sub _tempFiles{
	my $self=shift;
	$self->{'__tempFiles'}=shift // return $self->{'__tempFiles'};
}

sub _initialize{
	my $self=shift;

	#init values
	my %params = @_;

	$self->_baseDirectory($params{'baseDirectory'}) if defined $params{'baseDirectory'};
	
	#logging
	$self->logger(Log::Log4perl->get_logger());

	#set anonymous values
	$self->_tempFiles([]);
}

sub run{
	my $self=shift;
	
	my %settings = @_;
	$self->_queryFile($settings{'queryFile'}) // confess("queryFile required in MummerGPU");
	$self->_referenceFileArray([$settings{'referenceFile'}]) if defined $settings{'referenceFile'};
	$self->_mummerDirectory($settings{'mummerDirectory'}) // confess("mummerDirectory required in MummerGPU");
	$self->_baseDirectory($settings{'baseDirectory'}) // confess("baseDirectory required in MummerGPU");
	$self->_refFileSizeLimit($settings{'refFileSizeLimit'}) // $self->_refFileSizeLimit(100000000);
	$self->_numberOfCores($settings{'numberOfCores'}) // $self->_numberOfCores(1);
	$self->_gpu($settings{'gpu'}) // $self->_gpu(0);
	$self->_b($settings{'b'}) // $self->_b(200);
	$self->_c($settings{'c'}) // $self->_c(50);
	$self->_d($settings{'d'}) // $self->_d(0.12);
	$self->_l($settings{'l'}) // $self->_l(20);
	$self->_g($settings{'g'}) // $self->_g(100);
	$self->_p($settings{'p'}) // $self->_p($self->_baseDirectory . 'out');
	
	#check for defined reference files
	unless(defined $self->_referenceFileArray->[0]){
		$self->logger->fatal("MummerGPU required at least one reference file.\n Either run MummerGPU::mummersLittleHelper 
			or specify a reference file");
		exit(1);
	}

	$self->_systemLineBase($self->_mummerDirectory);
	if($self->_gpu == 1){
		$self->_systemLineBase($self->_systemLineBase . 'mummergpu');
	}
	else{
		$self->_systemLineBase($self->_systemLineBase. 'nucmer');
	}						
	$self->_systemLineBase($self->_systemLineBase . ' --maxmatch -b '. $self->_b . ' -c ' . $self->_c . ' -d ' . $self->_d . ' -g ' . $self->_g . ' -l ' . $self->_l);
	
	my $forkManager = Parallel::ForkManager->new($self->_numberOfCores);
	
	my $tempNum=0;	
	foreach my $refFile(@{$self->_referenceFileArray}){
		$tempNum++;	
		my $tempFileName = $self->_p . '_temp'.$tempNum;
		push @{$self->_tempFiles}, $tempFileName;

		$forkManager->start and next();				
			$self->_launchMummer(
				$tempFileName,
				$refFile,
				$self->_queryFile
			);				
		$forkManager->finish;			
	}	

	$forkManager->wait_all_children();
	$self->deltaFile($self->_p . '.delta');
	$self->logger->info("Delta file: " . $self->deltaFile);
	
	$self->_combineDeltaFiles();

	#remove the temp files
	if(defined $self->_referenceFileArray->[0]){
		$self->_removeSplitFiles();
	}
}

sub _removeSplitFiles{
	my $self=shift;

	foreach my $splitFile(@{$self->_referenceFileArray}){
		unlink $splitFile;
	}
}


sub _launchMummer{
	my $self=shift;
	
	#-p as first option,
	#refFileName as second option,
	#queryFileName as third option
	
	my $p=shift;
	my $ref=shift;
	my $query=shift;
	

	my $systemLine=$self->_systemLineBase . ' -p ' . $p . ' ' . $ref . ' ' . $query;
	$self->logger->info("Launching mummer with: $systemLine");
	system($systemLine);	
}

sub _combineDeltaFiles{
	my $self=shift;

	#combine the deltaFiles
	my $deltaFH=IO::File->new('>'.$self->deltaFile) or die "$!";
	$self->logger->info("Combining the delta files");		
	
	my $manny = FileInteraction::FileManipulation->new();
	$manny->outputFH($deltaFH);
	
	foreach my $file(@{$self->_tempFiles}){ 
		my $deltaFile = $file . '.delta';
		my $dFH = IO::File->new('<'. $deltaFile) or die "cannot open $deltaFile $!";
		while(my $line=$dFH->getline){
			$deltaFH->print($line);
		}
		$dFH->close();
		unlink $deltaFile;
	}	
	$deltaFH->close();
}

sub showCoords{
	my $self=shift;

	my %params = @_;
	my $systemLine = $self->_mummerDirectory . 'show-coords ' . $params{'deltaFile'} . ' -l -q -T > ' . $params{'coordsFile'};

	$self->logger->info("Launching show-coords with $systemLine");

	system($systemLine);
	$self->coordsFile($params{'coordsFile'});
}

sub mummersLittleHelper{
	my $self=shift;

	my %params = @_;

	if(defined $params{'bpPerFile'} && defined $params{'numberOfFiles'}){
		$self->logger->fatal("Only one of bpPerFile or numberOfFiles can be defined in mummersLittleHelper");
		exit(1);
	}

	#check that at least one is set by the user, or complain
	unless(defined $params{'bpPerFile'} || defined $params{'numberOfFiles'}){
		$self->logger->fatal("mummersLittleHelper requires either bpPerFile or numberOfFiles to be set");
		exit(1);
	}
	#if bpPerFile is set, create unlimited number of files at the threshold
	#if numberOfFiles is set, create fixed number of files with bp/file variable

	my $bpPerFile;
	if(defined $params{'bpPerFile'}){
		$bpPerFile = $params{'bpPerFile'};
	}
	else{
		my $multiFastaFileSize = (-s $params{'multiFastaFile'});
		$bpPerFile = int($multiFastaFileSize / $params{'numberOfFiles'})+1;
	}

	my $multiFastaFH = Bio::SeqIO->new(-file=>'<'.$params{'multiFastaFile'}, -format=>'fasta') or die "$!";

	my $toCreateNewFile=1;
	my $tempCounter=0;
	my $outputFH;
	my $bpSize=0;
	my @tempFastaFiles;

	while(my $seq = $multiFastaFH->next_seq()){
		if($toCreateNewFile == 1){
			unless($tempCounter ==0){
				$outputFH->close();
			}

			my $fileName = $self->_baseDirectory . 'multiFastaTemp_' . $tempCounter;
			$outputFH = Bio::SeqIO->new(-file=>'>'.$fileName , -format=>'fasta') or die "$!";
			$tempCounter++;
			$toCreateNewFile=0;
			$bpSize=0;
			push @tempFastaFiles,$fileName;
		}

		$bpSize +=$seq->length();
		$outputFH->write_seq($seq);

		if($bpSize >= $bpPerFile){
			$toCreateNewFile=1;
		}
	}

	$outputFH->close();
	$multiFastaFH->close();
	$self->_referenceFileArray(\@tempFastaFiles);
}

1;

