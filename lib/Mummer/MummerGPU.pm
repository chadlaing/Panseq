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
use Logging::Logger;
use File::Temp;
use Parallel::ForkManager;
use FileInteraction::FileManipulation;
our @ISA = qw/Logger/;

#object creation
sub new{
	my $class=shift;
    my $self = {};
    bless ($self, $class);
    $self->_createMummerFilesInitialize(@_);   
    return $self;
}

#class variables
sub _queryFile{
	my $self=shift;
	$self->{'_queryFile'}=shift // return $self->{'_queryFile'};
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


sub _createMummerFilesInitialize{
	my $self=shift;
	
	#inheritance
	$self->_loggerInitialize(@_);	
}

sub run{
	my $self=shift;
	
	my $settingsHashRef = shift;
	$self->_queryFile($settingsHashRef->{'queryFile'}) // confess("queryFile required in MummerGPU");
	$self->_referenceFile($settingsHashRef->{'referenceFile'}) // confess("referenceFile required in MummerGPU");
	$self->_mummerDirectory($settingsHashRef->{'mummerDirectory'}) // confess("mummerDirectory required in MummerGPU");
#	$self->deltaFile($settingsHashRef->{'deltaFile'}) // confess("deltaFile required in MummerGPU");
	$self->_baseDirectory($settingsHashRef->{'baseDirectory'}) // confess("baseDirectory required in MummerGPU");
	$self->_numberOfCores($settingsHashRef->{'numberOfCores'}) // $self->_numberOfCores(1);
	$self->_gpu($settingsHashRef->{'gpu'}) // $self->_gpu(0);
	$self->_b($settingsHashRef->{'b'}) // $self->_b(200);
	$self->_c($settingsHashRef->{'c'}) // $self->_c(50);
	$self->_d($settingsHashRef->{'d'}) // $self->_d(0.12);
	$self->_l($settingsHashRef->{'l'}) // $self->_l(20);
	$self->_g($settingsHashRef->{'g'}) // $self->_g(100);
	$self->_p($settingsHashRef->{'p'}) // $self->_p('out');
	
#	my $forkManager = Parallel::ForkManager->new($self->_numberOfCores);
#	{
#		local $/='>';
#		my $refFH = IO::File->new('<'. $self->_referenceFile) or die "$!";
#		
#		my $tempNum=0;		
#		while(<$refFH>){
#			my $fasta=$_;
#			next if $fasta eq '>';
#			
#			my $header;
#			my $seq;
#			if($fasta =~ /(^\N+)\R+/){
#				$header=$1;
#				$seq=$fasta;
#				$seq =~ s/\Q$header//;
#				$seq =~ s/\R//g;
#			}
#			else{
#				$self->logger->fatal("FATAL:\theader: $header");
#				$self->logger->fatal("FATAL:\tseq: $seq");
#				$self->logger->fatal("FATAL:\tfasta: $fasta");
#				exit(1);
#			}
#			
#			$forkManager->start and next();
#				#print out the currentfile as temp
#				my $tempRefFH = File::Temp->new();
#				my $tempQueryFH = File::Temp->new();
#				
#				#create ref temp
#				print $tempRefFH '>'. $header . "\n" . $seq . "\n";
#		
				#run the mummerGPU
				my $systemLine = $self->_mummerDirectory;
				
				if($self->_gpu == 1){
					$systemLine .='mummergpu';
			}
				else{
					$systemLine .='nucmer';
				}		
				
				$systemLine .= ' --maxmatch -b '. $self->_b . ' -c ' . $self->_c . ' -d ' . $self->_d . ' -g ' . $self->_g . ' -l ' . $self->_l . ' -p ' . $self->_baseDirectory . $self->_p . ' ' . $self->_referenceFile . ' ' . $self->_queryFile;
#				$systemLine .= ' -b '. $self->_b . ' -c ' . $self->_c . ' -d ' . $self->_d . ' -g ' . $self->_g . ' -l ' . $self->_l . ' -p ' . $self->_baseDirectory . $tempNum . '_temp'. ' ' . $tempRefFH->filename . ' ' . $self->_queryFile;
#			
#			$self->logger->info("INFO:\tsystemLine: $systemLine");			
			system($systemLine);
			$self->deltaFile($self->_baseDirectory . $self->_p . '.delta');
			$self->logger->info("INFO:\tDelta file: " . $self->deltaFile);
#			$forkManager->finish;			
#			
#		}
#		continue{ 
#			$tempNum++;
#		}
#		$refFH->close();		
#	}
#	
#	$forkManager->wait_all_children();
#	$self->_combineDeltaFiles();
}


sub _combineDeltaFiles{
	my $self=shift;
	
	#combine the deltaFiles
	my $deltaFH=IO::File->new('>'.$self->deltaFile) or die "$!";
	$self->logger->info("INFO:\tCombining the delta files");		
	
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

