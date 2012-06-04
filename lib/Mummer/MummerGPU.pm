#!/usr/bin/perl

#takes in the queryFile (all query files in one file)
#and the referenceFile (all reference files in one file)
#mummer needs the single query file with all seqs vs each ref sequence in a separate file

package MummerGPU;

use strict;
use warnings;
use FindBin::libs;
use Carp;
use IO::File;
use File::Temp;
use Parallel::ForkManager;
use FileInteraction::FileManipulation;
use Log::Log4perl;
use Bio::SeqIO;
use File::Copy;

#object creation
sub new{
	my $class=shift;
    my $self = {};
    bless ($self, $class);
    $self->_mummerGPUInitialize(@_);   
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

sub _referenceFile{
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


sub _mummerGPUInitialize{
	my $self=shift;
	
	#logging
	$self->logger(Log::Log4perl->get_logger());
}

sub run{
	my $self=shift;
	
	my $settingsHashRef = shift;
	$self->_queryFile($settingsHashRef->{'queryFile'}) // confess("queryFile required in MummerGPU");
	$self->_referenceFile($settingsHashRef->{'referenceFile'}) // confess("referenceFile required in MummerGPU");
	$self->_mummerDirectory($settingsHashRef->{'mummerDirectory'}) // confess("mummerDirectory required in MummerGPU");
	$self->_baseDirectory($settingsHashRef->{'baseDirectory'}) // confess("baseDirectory required in MummerGPU");
	$self->_refFileSizeLimit($settingsHashRef->{'refFileSizeLimit'}) // $self->_refFileSizeLimit(100000000);
	$self->_numberOfCores($settingsHashRef->{'numberOfCores'}) // $self->_numberOfCores(1);
	$self->_gpu($settingsHashRef->{'gpu'}) // $self->_gpu(0);
	$self->_b($settingsHashRef->{'b'}) // $self->_b(200);
	$self->_c($settingsHashRef->{'c'}) // $self->_c(50);
	$self->_d($settingsHashRef->{'d'}) // $self->_d(0.12);
	$self->_l($settingsHashRef->{'l'}) // $self->_l(20);
	$self->_g($settingsHashRef->{'g'}) // $self->_g(100);
	$self->_p($settingsHashRef->{'p'}) // $self->_p('out');
	
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
	my $tempRefString='';
	my $refFileSize=0;		
	my $clearFlag=0;	

	my $refFastaFH = Bio::SeqIO->new(
                            -file   => $self->_referenceFile,
                            -format => 'fasta'
                    );   			
	
	while(my $fasta = $refFastaFH->next_seq()){			
			$clearFlag=0;
			
			$tempRefString .='>' . $fasta->id . $fasta->desc . "\n" . $fasta->seq . "\n";
			$refFileSize += $fasta->length;				
			
			if($refFileSize < $self->_refFileSizeLimit){
				next;
			}
			else{
				$self->logger->info("Temp ref file size (bp): $refFileSize");
				$refFileSize=0;
			}
						
			$forkManager->start and next();
				#print out the currentfile as temp
				my $tempRefFH = File::Temp->new();
				#my $tempQueryFH=File::Temp->new();
				#copy($self->_queryFile, $tempQueryFH->filename());
				
				$self->logger->debug("Temp ref file name: " . $tempRefFH->filename);
				#create ref temp
				print $tempRefFH $tempRefString;
		
				$self->_launchMummer(
					$self->_baseDirectory . $tempNum . '_temp',
					$tempRefFH->filename,
					#$tempQueryFH->filename()
					$self->_queryFile
				);				
			$forkManager->finish;			
			
		}
		continue{ 
			$tempNum++;
			if($clearFlag==1){
				$tempRefString='';
			}
		}
	$refFastaFH->close();		

	#run mummer on any stored values
	if($refFileSize > 0){
		my $tempRefFH = File::Temp->new();
		print $tempRefFH $tempRefString;
		$self->logger->info("Temp ref file size (bp): $refFileSize");
		
		$self->_launchMummer(
			$self->_baseDirectory . $tempNum . '_temp',
			$tempRefFH->filename,
			$self->_queryFile
		);	
		$tempRefFH->close();
	}	
	$forkManager->wait_all_children();
	$self->deltaFile($self->_baseDirectory . $self->_p . '.delta');
	$self->logger->info("Delta file: " . $self->deltaFile);
	
	$self->_combineDeltaFiles();
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
	
	my $manny = FileManipulation->new();
	$manny->outputFilehandle($deltaFH);
	my $fileNamesRef = $manny->getFileNamesFromDirectory($self->_baseDirectory);
	
	foreach my $file(@{$fileNamesRef}){ 

		unless($file =~ m/_temp\.delta/){
			next;
		}
		
		my $dFH = IO::File->new('<'.$file) or die "$!";
		while(my $line=$dFH->getline){
			print $deltaFH $line;
		}
		$dFH->close();
		unlink $file;
	}	
	$deltaFH->close();
}

1;

