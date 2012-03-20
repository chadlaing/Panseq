#!/usr/bin/perl
package SequenceRetriever;

use strict;
use warnings;
use Carp;
use FindBin::libs;
use IO::File;
use Bio::DB::Fasta;
use FileInteraction::FlexiblePrinter;
our @ISA = qw/FlexiblePrinter/;

#object creation
use Object::Tiny::RW qw{
	databaseName
	formattedFileName
};

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_sequenceRetrieverInitialize(@_);
    return $self;
}


##methods
sub _sequenceRetrieverInitialize{
	my($self)=shift;
    my $file=shift;
   
    #inheritance
    $self->_flexiblePrinterInitialize(@_);
    
    if($file){    	
    	$self->createDatabase($file);
    }
    else{
    	print STDERR "no database file specified!\n";
    	exit(1);
    }
}

sub createDatabase{
	my($self)=shift;
	
	if(@_){
		my $file=shift;
		$self->databaseName(Bio::DB::Fasta->new($self->createProperlyFormattedFastaFile($file)));
	}
	else{
		print STDERR "databaseFile not defined!\n";
		exit(1);
	}
}

sub createProperlyFormattedFastaFile{
	my($self)=shift;
	
	if(@_){
		my $file=shift;

		open (my $originalfile,'<',$file) or die "$!";
		my $formattedFileName = $file . '_dbTEMP';
		#store
		$self->formattedFileName($formattedFileName);
		
		my $oFile = IO::File->new('>' . $formattedFileName);
		local $/='>';
		
		while(my $segment = <$originalfile>){
			next if $segment eq '>';
			
			my $header;
			my $sequence;
			
			if($segment =~ /(^.+)/){
				$header=$1;
				$segment =~ s/\Q$header//;
				$segment =~ s/\W//g;
				$sequence = $segment;
				undef $segment;
			}
			else{
				print STDERR "not a fasta fragment in createProperlyFormattedFastaFile!\n";
				print STDERR "segment:\n$segment\n";
				exit(1);
			}
			
			print $oFile ('>' . $header . "\n");
			my $maxChars=65534;
			
			for(my $i=0; $i<length($sequence); $i+=$maxChars){
				my $charsLeft = length($sequence)-$i +1;
				my $extractionLength;
				
				if($charsLeft >= $maxChars){
					$extractionLength=$maxChars;
				}
				else{
					$extractionLength = $charsLeft;
				}
				
				my $subseq = substr($sequence,$i,$extractionLength);
				print $oFile ($subseq . "\n");
			}
		}
		close $originalfile;
		$oFile->close();
		return $formattedFileName;
	}
	else{
		print STDERR "nothing sent to createProperlyFormattedFastaFile\n";
		exit(1);
	}
}


sub extractAndPrintRegionsFromHash{
	my($self)=shift;
	
	my $paramsRef=shift;
	
	my $hashRef = $paramsRef->{'novelRegionHashRef'} // confess ('Hash reference required in extractAndPrintRegionsFromHash');
	my $cutoffSize = $paramsRef->{'cutoffSize'} // confess ('cutoffSize required in extractAndPrintRegionsFromHash');
	$self->outputFilehandle($paramsRef->{'novelOutputFH'}) if defined $paramsRef->{'novelOutputFH'};
	

	#expects hashvalues of $hash->id->,1..5,8..134,678..45999 etc. 
	#same format as generated from NovelRegionFinder
		
	foreach my $id(keys %{$hashRef}){
		my $coordString = $hashRef->{$id};
			
		my $novelCounter=0;
		while($coordString =~ /\,(\d+)\.\.(\d+)/gc){				
			my $start =$1;
			my $end =$2;
			my $length=$end-$start+1;
				
			if(defined $cutoffSize){
				next unless $length >= $cutoffSize;
			}
				
			$novelCounter++;
				
			$self->printOut('>' . $id . '|NovelRegion=' . $novelCounter . '|Start=' . $start . '|End=' . $end . '|Length=' . $length . "\n" . $self->extractRegion($id,$start,$end) . "\n");	
		}# end while
	}# end foreach	
}

sub extractRegion{
	my($self)=shift;
	
	if(@_){
		my $seqID=shift;
		my $startBp;
		my $endBp;
		
		if(scalar(@_)==2){
			$startBp=shift;
			$endBp=shift;
		}

		my $seqObj = $self->databaseName->get_Seq_by_id($seqID) // confess ("Cannot locate $seqID");
		
		if(defined $startBp && defined $endBp){
			return $seqObj->subseq($startBp=>$endBp);
		}
		else{
			return $seqObj->seq;
		}					
	}
	else{
		print STDERR "extract region is missing parameters!\n";
		exit(1);
	}
}

1;
