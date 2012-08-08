#!/usr/bin/perl

package Blast::BlastParallel;

use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin";
use IO::File;
use Parallel::ForkManager;
use FileInteraction::Fasta::FastaFileSplitter;
use Log::Log4perl;

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_initialize(@_);
    return $self;
}

#class variables
sub arrayOfXMLFiles{
	my $self=shift;
	$self->{'_arrayOfXMLFiles'} = shift // return $self->{'_arrayOfXMLFiles'}
}

sub combinedXMLFileName{
	my $self=shift;
	$self->{'_combinedXMLFileName'} = shift // return $self->{'_combinedXMLFileName'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

sub _initialize{
	my($self)=shift;
	
	#logging
	$self->logger(Log::Log4perl->get_logger());
	
	#anonymous array
	$self->arrayOfXMLFiles([]);
}
sub runBlastParallel{
	my($self)=shift;
	
	my $paramsRef=shift;
	
	my $blastIO=$paramsRef->{'blastIO'} // confess ('blastIO objectRequired in runBlastParallel');
	my $inputFile=$paramsRef->{'inputLociFile'} // confess ('inputLociFile required in runBlastParallel');
	my $numberOfCores=$paramsRef->{'numberOfCores'}  // confess ('numberOfCores required in runBlastParallel');
	
	my $splitter = FastaFileSplitter->new();
	$splitter->splitFastaFile($inputFile,$numberOfCores);
		
	my $forker = Parallel::ForkManager->new($numberOfCores);

	foreach my $blastFile(@{$splitter->arrayOfSplitFiles}){
		my $blastOutputFile= $blastFile . '_blastoutput.xml';
			
		$self->logger->debug("DEBUG:\tCreating temp xml file: $blastOutputFile");
			
		push @{$self->arrayOfXMLFiles}, $blastOutputFile;			
		$forker->start and next;
			$blastIO->runBlastn(
				'query'=>$blastFile,
				'out'=>$blastOutputFile
			);	
		unlink $blastFile;
		#end run blast
		$forker->finish;				
	}
	$forker->wait_all_children();	
}

sub combineXMLFiles{
	my($self)=shift;
	
	if(@_){
		my $directory=shift;
		
		#create master XML file
		my $combiner = FileManipulation->new();
		my $finalBlastXMLfileName = $directory . 'masterBlast.xml';
		my $finalBlastXML = IO::File->new('>'. $finalBlastXMLfileName) or die "$!";
		$combiner->outputFilehandle($finalBlastXML);		
		$combiner->vanillaCombineFiles($self->arrayOfXMLFiles,1); #true for unlink
		$finalBlastXML->close();		
		$self->combinedXMLFileName($finalBlastXMLfileName);
	}
	else{
		print STDERR "no directory specified in combineTempFiles\n";
		exit(1);
	}	
}

1;
	