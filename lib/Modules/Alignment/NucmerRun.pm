#!/usr/bin/env perl

package Modules::Alignment::NucmerRun;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Log::Log4perl;
use Carp;

#object creation
sub new{
	my $class=shift;
    my $self = {};
    bless ($self, $class);
    $self->_initialize(@_);   
    return $self;
}

sub _initialize{
	my $self=shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::NucmerRun\n");

	#init values
	my %params = @_;

	foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			$self->logger->logconfess("$key is not a valid parameter in Modules::Alignment::NucmerRun");
		}	
	}
}

#class variables
sub logFile{
	my $self=shift;
	$self->{'_mummerLog'} = shift // return $self->{'_mummerLog'};
}

sub mummerDirectory{
	my $self=shift;
	$self->{'_mummerDirectory'} = shift // return $self->{'_mummerDirectory'};
}

sub queryFile{
	my $self=shift;
	$self->{'_queryFile'}=shift // return $self->{'_queryFile'};
}

sub referenceFile{
	my $self=shift;
	$self->{'_referenceFile'}=shift // return $self->{'_referenceFile'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub numberOfCores{
	my $self=shift;
	$self->{'_numberOfCores'}=shift // return $self->{'_numberOfCores'};
}

sub deltaFile{
	my $self=shift;
	$self->{'_deltaFile'}=shift // return $self->{'_deltaFile'};
}

sub referenceFileArray{
	my $self=shift;
	$self->{'_referenceFile'}=shift // return $self->{'_referenceFile'};
}

sub refFileSizeLimit{
	my $self=shift;
	$self->{'_refFileSizeLimit'}=shift // return $self->{'_refFileSizeLimit'};	
}

sub systemLineBase{
	my $self=shift;
	$self->{'_systemLineBase'}=shift // return $self->{'_systemLineBase'};
}


sub l{
	my $self=shift;
	$self->{'_l'}=shift // return $self->{'_l'};
}

sub g{
	my $self=shift;
	$self->{'_g'}=shift // return $self->{'_g'};
}

sub d{
	my $self=shift;
	$self->{'_d'}=shift // return $self->{'_d'};
}

sub c{
	my $self=shift;
	$self->{'_c'}=shift // return $self->{'_c'};
}

sub b{
	my $self=shift;
	$self->{'_b'}=shift // return $self->{'_b'};
}

sub p{
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



sub run{
	my $self=shift;
	
	#check for required files
	unless(defined $self->referenceFile){
		$self->logger->logconfess("Reference file required in run\n");
	}

	unless(defined $self->queryFile){
		$self->logger->logconfess("Query file required in run\n");
	}

	unless(defined $self->p){
		$self->logger->logconfess("Delta file prefix required in run\n");
	}

	unless(defined $self->mummerDirectory){
		$self->logger->logconfess("Mummer directory required in run");
	}

	#run mummer
	my $mummerLine = $self->_createMummerLine();
	$self->logger->info("Running nucmer comparison with the command: $mummerLine");
	system($mummerLine);	

	#add a delta-filter to limit the nucmer matches to the percentIdentity cutoff

	#run show-coords on delta file
	$self->_showCoords();
}

sub _createMummerLine{
	my $self=shift;

	my $mummerLine = $self->mummerDirectory . 'nucmer --maxmatch';

	$mummerLine .= ' -b ' . $self->b if defined $self->b;
	$mummerLine .= ' -c ' . $self->c if defined $self->c;
	$mummerLine .= ' -d ' . $self->d if defined $self->d;
	$mummerLine .= ' -g ' . $self->g if defined $self->g;
	$mummerLine .= ' -l ' . $self->l if defined $self->l;
	$mummerLine .= ' -b ' . $self->b if defined $self->b;
	$mummerLine .= ' -p ' . $self->p if defined $self->p;
	$mummerLine .= ' ' . $self->referenceFile;
	$mummerLine .= ' ' . $self->queryFile;

	return $mummerLine;
}


sub _showCoords{
	my $self=shift;

	#check for requirements
	unless(defined $self->coordsFile){
		$self->logger->logconfess("Coords file name required in _showCoords");
	}

	my $coordsLine = $self->mummerDirectory . 'show-coords ' . $self->p . '.delta' . ' -l -q -T > ' . $self->coordsFile;

	$self->logger->info("Launching show-coords with $coordsLine");

	system($coordsLine);
}

# sub _removeSplitFiles{
# 	my $self=shift;

# 	foreach my $splitFile(@{$self->_referenceFileArray}){
# 		unlink $splitFile;
# 	}
# }


# sub _combineDeltaFiles{
# 	my $self=shift;

# 	#combine the deltaFiles
# 	my $deltaFH=IO::File->new('>'.$self->deltaFile) or die "$!";
# 	$self->logger->info("Combining the delta files");		
	
# 	my $manny = FileInteraction::FileManipulation->new();
# 	$manny->outputFH($deltaFH);
	
# 	foreach my $file(@{$self->_tempFiles}){ 
# 		my $deltaFile = $file . '.delta';
# 		my $dFH = IO::File->new('<'. $deltaFile) or die "cannot open $deltaFile $!";
# 		while(my $line=$dFH->getline){
# 			$deltaFH->print($line);
# 		}
# 		$dFH->close();
# 		unlink $deltaFile;
# 	}	
# 	$deltaFH->close();
# }

# sub mummersLittleHelper{
# 	my $self=shift;

# 	my %params = @_;

# 	if(defined $params{'bpPerFile'} && defined $params{'numberOfFiles'}){
# 		$self->logger->fatal("Only one of bpPerFile or numberOfFiles can be defined in mummersLittleHelper");
# 		exit(1);
# 	}

# 	#check that at least one is set by the user, or complain
# 	unless(defined $params{'bpPerFile'} || defined $params{'numberOfFiles'}){
# 		$self->logger->fatal("mummersLittleHelper requires either bpPerFile or numberOfFiles to be set");
# 		exit(1);
# 	}
# 	#if bpPerFile is set, create unlimited number of files at the threshold
# 	#if numberOfFiles is set, create fixed number of files with bp/file variable

# 	my $bpPerFile;
# 	if(defined $params{'bpPerFile'}){
# 		$bpPerFile = $params{'bpPerFile'};
# 	}
# 	else{
# 		my $multiFastaFileSize = (-s $params{'multiFastaFile'});
# 		$bpPerFile = int($multiFastaFileSize / $params{'numberOfFiles'})+1;
# 	}

# 	my $multiFastaFH = Bio::SeqIO->new(-file=>'<'.$params{'multiFastaFile'}, -format=>'fasta') or die "$!";

# 	my $toCreateNewFile=1;
# 	my $tempCounter=0;
# 	my $outputFH;
# 	my $bpSize=0;
# 	my @tempFastaFiles;

# 	while(my $seq = $multiFastaFH->next_seq()){
# 		if($toCreateNewFile == 1){
# 			unless($tempCounter ==0){
# 				$outputFH->close();
# 			}

# 			my $fileName = $self->_baseDirectory . 'multiFastaTemp_' . $tempCounter;
# 			$outputFH = Bio::SeqIO->new(-file=>'>'.$fileName , -format=>'fasta') or die "$!";
# 			$tempCounter++;
# 			$toCreateNewFile=0;
# 			$bpSize=0;
# 			push @tempFastaFiles,$fileName;
# 		}

# 		$bpSize +=$seq->length();
# 		$outputFH->write_seq($seq);

# 		if($bpSize >= $bpPerFile){
# 			$toCreateNewFile=1;
# 		}
# 	}

# 	$outputFH->close();
# 	$multiFastaFH->close();
# 	$self->_referenceFileArray(\@tempFastaFiles);
# }

1;