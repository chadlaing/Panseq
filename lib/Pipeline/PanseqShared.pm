#!/usr/bin/perl

package PanseqShared;

use strict;
use warnings;
use diagnostics;
use FindBin::libs;
use Carp;
use File::Path qw{make_path remove_tree};
use File::Copy;
use Logging::Logger;
use FileInteraction::LinePackage;
our @ISA = qw/Logger/;

#object creation
sub new{
	my($class)=shift;
    my $self = {};
    bless ($self, $class);
    $self->_panseqSharedInitialize(@_);
    return $self;
}

#class variables
sub queryDirectory{
	my $self=shift;
	$self->{'_PanseqShared_queryDirectory'}=shift // return $self->{'_PanseqShared_queryDirectory'};
}

sub referenceDirectory{
	my $self=shift;
	$self->{'_PanseqShared_referenceDirectory'}=shift // return $self->{'_PanseqShared_referenceDirectory'};
}




sub _baseDirectory{
	my $self=shift;
	$self->{'_PanseqShared_baseDirectory'}=shift // return $self->{'_PanseqShared_baseDirectory'};
}

sub _numberOfCores{
	my $self=shift;
	$self->{'_PanseqShared_numberOfCores'}=shift // return $self->{'_PanseqShared_numberOfCores'};
}

sub mummerDirectory{
	my $self=shift;
	$self->{'_PanseqShared_mummerDirectory'}=shift // return $self->{'_PanseqShared_mummerDirectory'};
}

sub minimumNovelRegionSize{
	my $self=shift;
	$self->{'_PanseqShared_minimumNovelRegionSize'}=shift // return $self->{'_PanseqShared_minimumNovelRegionSize'};
}


sub createGraphic{
	my $self=shift;
	$self->{'_PanseqShared_createGraphic'}=shift // return $self->{'_PanseqShared_createGraphic'};
}

sub combinedQueryFile{
	my $self=shift;
	$self->{'_PanseqShared_combinedQueryFile'}=shift // return $self->{'_PanseqShared_combinedQueryFile'};
}

sub queryNameObjectHash{
	my $self=shift;
	$self->{'_PanseqShared_queryNameObjectHash'}=shift // return $self->{'_PanseqShared_queryNameObjectHash'};
}

sub _novel_nucB{
	my $self=shift;
	$self->{'_PanseqShared_novel_nucB'}=shift // return $self->{'_PanseqShared_novel_nucB'};
}

sub _novel_nucC{
	my $self=shift;
	$self->{'_PanseqShared_novel_nucC'}=shift // return $self->{'_PanseqShared_novel_nucC'};
}

sub _novel_nucD{
	my $self=shift;
	$self->{'_PanseqShared_novel_nucD'}=shift // return $self->{'_PanseqShared_novel_nucD'};
}

sub _novel_nucG{
	my $self=shift;
	$self->{'_PanseqShared_novel_nucG'}=shift // return $self->{'_PanseqShared_novel_nucG'};
}

sub _novel_nucL{
	my $self=shift;
	$self->{'_PanseqShared_novel_nucL'}=shift // return $self->{'_PanseqShared_novel_nucL'};
}

#methods
sub _panseqSharedInitialize{
	my $self=shift;
	$self->_loggerInitialize(@_);
}

sub getSettingsFromConfigurationFile{
	my($self)=shift;
	
	if(@_){
		my $fileName=shift;
		my $inFile = IO::File->new('<' . $fileName) or die "\nCannot open $fileName\n $!";
		my %settingsHash;
		
		while(my $line = $inFile->getline){
			next unless $line =~ m/\t/;
			my $la = LinePackage->new($line);
			
			my $setting = $la->lineArray->[0] // 'no_value';
			my $value = $la->lineArray->[1] // 'no_value';
			$settingsHash{$setting}=$value;	
		}		
		$inFile->close();
		
		$self->validateUniversalSettings(\%settingsHash);
		return \%settingsHash;
	}
	else{
		confess "no file specified in getSettingsFromConfigurationFile\n";
	}
}

sub missingParam{
	my($self)=shift;
	
	print STDERR "Missing parameter $_[0], which is required for a NovelRegionFinder run!\n";
	exit(1);
}

sub nucDCheck{
	my($self)=shift;
	
	if(@_){
		my $type=shift;
		
		if($type =~ /^0\.\d+$/){
			if($type > 0){
				return $type;
			}
			else{
				confess "nucD must be greater than 0 in configuration file!\n";
			}			
		}
		else{
			confess "$type is an invalid entry for nucD in configuration file!\n",
				"Valid types are real number values 0 < X < 1\n";
		}
	}
	else{
		confess "Nothing sent to nucDCheck!\n";
	}
}

sub validateUniversalSettings{
	my($self)=shift;
	
	if(@_){		
		my $settingsHashRef=shift;		
		
		foreach my $setting(keys %{$settingsHashRef}){
			my $value = $settingsHashRef->{$setting};
			
			#directories
			$self->queryDirectory($self->isADirectory($value)) if $setting eq 'queryDirectory';
			$self->_baseDirectory($self->isADirectory($value)) if $setting eq 'baseDirectory';
			$self->mummerDirectory($self->isADirectory($value)) if $setting eq 'mummerDirectory';
			$self->referenceDirectory($self->isADirectory($value)) if $setting eq 'referenceDirectory';

			#integers
			$self->_numberOfCores($self->isAnInt($value)) if $setting eq 'numberOfCores';	
			$self->minimumNovelRegionSize($self->isAnInt($value)) if $setting eq 'minimumNovelRegionSize';
			
			#yes or no
			$self->createGraphic($self->yesOrNoCheck($value)) if $setting eq 'createGraphic';		
		}
		
		#defaults and requirements
		$self->missingParam('queryDirectory') unless (defined $self->queryDirectory);
		$self->missingParam('baseDirectory') unless defined $self->_baseDirectory;
		$self->missingParam('mummerDirectory') unless defined $self->mummerDirectory;
		$self->_numberOfCores(1) unless defined $self->_numberOfCores;
		$self->missingParam('minimumNovelRegionSize') unless defined $self->minimumNovelRegionSize;
		$self->createGraphic('no') unless defined $self->createGraphic;
	}
}


sub percentIdentityCheck{
	my($self)=shift;
	
	if(@_){
		my $number = shift;
		
		if(($self->isAnInt($number)) && ($number <= 100)){
			return $number;
		}
		else{
			confess "$number is an invalid entry for percentIdentityCutoff.\n",
				"Values must be integers 0 <= X <= 100";
		}
	}
	else{
		confess "Nothing sent to percentIdentityCheck!\n";
	}
}


sub yesOrNoCheck{
	my($self)=shift;
	
	if(@_){
		my $type=shift;
		
		if(($type eq 'yes') || ($type eq 'no')){
			return $type;
		}
		else{
			confess "$type is an invalid type!\n",
				"Valid types are yes or no.\n";
		}
	}
	else{
		confess "nothing sent to createGraphicCheck!\n";
	}
}


sub isAnInt{
	my($self)=shift;
	
	if(@_){
		my $value=shift;
		
		if($value =~ /^\d+$/){
			if($value > 0){
				return $value;
			}
			else{
				confess "Integer value required to be >0 in configuration file!\n";
			}
		}
		else{
			confess "$value is not an integer and should be in configuration file!\n";
		}
	}
	else{
		confess "nothing sent to isAnInt!\n";
	}
}

sub isADirectory{
	my($self)=shift;
	
	if(@_){
		my $dirName = shift;
		
		if($dirName =~ /^\//){
			unless($dirName =~ /\/$/){
				#add trailing slash if absent
				$dirName .= '/';
			}
			return $dirName;
		}
		else{
			confess "directory names must start with '/' in configuration file!\n";
		}
	}
	else{
		confess "nothing sent to isADirectory!\n";
	}
}

sub createDirectories{
	my($self)=shift;
	
	if(defined $self->_baseDirectory){
		
		$self->logger->info("INFO:\tCreating directories " . $self->_baseDirectory);
		#uses File::Path	
		remove_tree($self->_baseDirectory);
		make_path($self->_baseDirectory);
	}
	else{
		print STDERR "baseDirectory undefined!\n";
		exit(1);
	}
}

sub getQueryNamesAndCombineAllInputFiles{
	my($self)=shift;
	
	if(1){
		#get list of the query files
		#we need to create a queryFiles and referenceFiles file
		my $fileManipulator = FileManipulation->new();
		my $queryFileName = $self->_baseDirectory . 'queryFilesCombined.fasta';
		my $queryOutputFH = IO::File->new('>' .$queryFileName) or die "$!";
		$fileManipulator->outputFilehandle($queryOutputFH); #set filehandle
		
		$self->logger->info("INFO:\tCreating combined query files in $queryFileName");

		#store for posterity
		$self->combinedQueryFile($queryFileName);
		$fileManipulator->combineFilesIntoSingleFile($fileManipulator->getFileNamesFromDirectory($self->queryDirectory),1); #second argument cleans file
		$queryOutputFH->close();
		
		
		if($self->can('novelRegionFinderMode')){
			my $referenceFileName = $self->_baseDirectory . 'referenceFilesCombined.fasta';
		
			$self->logger->info("INFO:\tCreating combined reference files in $referenceFileName");
			
			$self->combinedReferenceFile($referenceFileName); #save for posterity
			my $refFH = IO::File->new('>'. $referenceFileName) or die "$!";
			$fileManipulator->outputFilehandle($refFH);			 
			$fileManipulator->combineFilesIntoSingleFile($fileManipulator->getFileNamesFromDirectory($self->referenceDirectory),1);
			$refFH->close();	
			
			#for a unique novelRegionFinder run, need an all vs. all comparison
			$self->logger->debug("DEBUG:\tnovelRegionFinderMode: " . $self->novelRegionFinderMode);
			
			if($self->novelRegionFinderMode eq 'common_to_all' || $self->novelRegionFinderMode eq 'unique'){
				my $outRefFH = IO::File->new('>>'. $self->combinedReferenceFile) or die "$!";
				$fileManipulator->outputFilehandle($outRefFH);
				$fileManipulator->vanillaCombineFiles([$self->combinedQueryFile]);
				$outRefFH->close();
				$self->logger->info("INFO:\tCreating all strains file as: " . $self->combinedReferenceFile);				
			}
			
			#unique requires that all vs. all be performed
			if($self->novelRegionFinderMode eq 'unique'){
				#with File::Copy
				copy($self->combinedReferenceFile, $self->combinedQueryFile) or die "Could not complete File::Copy $!";		
			}				
		}	
		
		$self->queryNameObjectHash($self->getSequenceNamesAsHashRef($fileManipulator->getFastaHeadersFromFile($self->combinedQueryFile)));
			
		$self->logger->info("INFO:\tGathered query names from $queryFileName");
	}
	else{
		print STDERR "The universe has reversed polarity\n";
		exit(1);
	}
}

sub getSequenceNamesAsHashRef{
	my($self)=shift;
	
	if(@_){
		my $fastaHeadersArrayRef=shift;
		
		#get query names		
		my %seqNames;
		foreach my $header(@{$fastaHeadersArrayRef}){			
			my $seqName = SequenceName->new($header);
			
			if(defined $seqNames{$seqName->name}){
				my $tempSeq = $seqNames{$seqName->name};
				$tempSeq->addToArrayOfHeaders($header);
				$seqNames{$seqName->name}=$tempSeq;
			}
			else{
				$seqNames{$seqName->name}=$seqName;
			}
			$self->logger->debug("DEBUG:\tAdding" .$seqName->name . "to QNOH");
		}#end foreach
		return \%seqNames;		
	}# end if
	else{
		print STDERR "nothing sent to getQuerySequenceNamesAsHashRef\n";
		exit(1);
	}
}


1;