#!/usr/bin/perl
package Visual::SvgDrawer;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../";
use FileInteraction::Fasta::SequenceName;

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;
	$self->_initialize(@_);
	return $self;
}

sub _initialize{
	my $self = shift;

	#init anonymous data structures
	$self->_sequenceSizes({});
}

sub _sequenceSizes{
	my $self=shift;
	$self->{'__sequenceSizes'} = shift // return $self->{'__sequenceSizes'};
}

sub _sequenceData{
	my $self=shift;
	$self->{'__sequenceData'} = shift // return $self->{'__sequenceData'};
}

sub _biggestSequence{
	my $self=shift;
	$self->{'__biggestSequence'} = shift // return $self->{'__biggestSequence'};
}

sub _panGenomeSize{
	my $self=shift;
	$self->{'__panGenomeSize'} = shift // return $self->{'__panGenomeSize'};
}

sub _allSequences{
	my $self=shift;
	$self->{'__allSequences'} = shift // return $self->{'__allSequences'};
}

sub _gapSize{
	my $self=shift;
	$self->{'__gapSize'} = shift // return $self->{'__gapSize'};
}

sub _gapCharacter{
	my $self=shift;
	$self->{'__gapCharacter'} = shift // return $self->{'__gapCharacter'};
}

sub _printCharacter{
	my $self=shift;
	$self->{'__printCharacter'} = shift // return $self->{'__printCharacter'};
}

sub _maxRefPosition{
	my $self=shift;
	$self->{'__maxRefPosition'} = shift // return $self->{'__maxRefPosition'};
}

sub _inputFile{
	my $self=shift;
	$self->{'__inputFile'} = shift // return $self->{'__inputFile'};
}

sub _outputFile{
	my $self=shift;
	$self->{'__outputFile'} = shift // return $self->{'__outputFile'};
}

sub _fastaQueryFile{
	my $self=shift;
	$self->{'__fastaQueryFile'} = shift // return $self->{'__fastaQueryFile'};
}

sub _queryLength{
	my $self=shift;
	$self->{'__queryLength'} = shift // return $self->{'__queryLength'};
}

sub _refLength{
	my $self=shift;
	$self->{'__refLength'} = shift // return $self->{'__refLength'};
}

sub _currentFileLine{
	my $self=shift;
	$self->{'__currentFileLine'} = shift // return $self->{'__currentFileLine'};
}

sub _queryStart{
	my $self=shift;
	$self->{'__queryStart'} = shift // return $self->{'__queryStart'};
}

sub _queryEnd{
	my $self=shift;
	$self->{'__queryEnd'} = shift // return $self->{'__queryEnd'};
}

sub _refStart{
	my $self=shift;
	$self->{'__refStart'} = shift // return $self->{'__refStart'};
}

sub _refEnd{
	my $self=shift;
	$self->{'__refEnd'} = shift // return $self->{'__refEnd'};
}

sub _refID{
	my $self=shift;
	$self->{'__refID'} = shift // return $self->{'__refID'};
}

sub _previousRefID{
	my $self=shift;
	$self->{'__previousRefID'} = shift // return $self->{'__previousRefID'};
}

sub _previousQueryID{
	my $self=shift;
	$self->{'__previousQueryID'} = shift // return $self->{'__previousQueryID'};
}

sub _queryID{
	my $self=shift;
	$self->{'__queryID'} = shift // return $self->{'__queryID'};
}

sub _percentID{
	my $self=shift;
	$self->{'__percentID'} = shift // return $self->{'__percentID'};
}

sub _minimumMatchLength{
	my $self=shift;
	$self->{'__minimumMatchLength'} = shift // return $self->{'__minimumMatchLength'};
}

sub _refSequenceLength{
	my $self=shift;
	$self->{'__refSequenceLength'} = shift // return $self->{'__refSequenceLength'};
}

sub _querySequenceLength{
	my $self=shift;
	$self->{'__querySequenceLength'} = shift // return $self->{'__querySequenceLength'};
}

sub _newID{
	my $self=shift;
	$self->{'__newID'} = shift // return $self->{'__newID'};
}

sub _imageWidth{
	my $self=shift;
	$self->{'__imageWidth'} = shift // return $self->{'__imageWidth'};
}

sub _imageHeight{
	my $self=shift;
	$self->{'__imageHeight'} = shift // return $self->{'__imageHeight'};
}

sub _graphicWidth{
	my $self=shift;
	$self->{'__graphicWidth'} = shift // return $self->{'__graphicWidth'};
}

sub _valueY{
	my $self=shift;
	$self->{'__valueY'} = shift // return $self->{'__valueY'};
}

sub _valueX{
	my $self=shift;
	$self->{'__valueX'} = shift // return $self->{'__valueX'};
}

sub _incrementY{
	my $self=shift;
	$self->{'__incrementY'} = shift // return $self->{'__incrementY'};
}

sub _bpPerUnit{
	my $self=shift;
	$self->{'__bpPerUnit'} = shift // return $self->{'__bpPerUnit'};
}

sub _imageUnits{
	my $self=shift;
	$self->{'__imageUnits'} = shift // return $self->{'__imageUnits'};
}

sub _rectWidth{
	my $self=shift;
	$self->{'__rectWidth'} = shift // return $self->{'__rectWidth'};
}

sub _rectHeight{
	my $self=shift;
	$self->{'__rectHeight'} = shift // return $self->{'__rectHeight'};
}

sub _rectColor{
	my $self=shift;
	$self->{'__rectColor'} = shift // return $self->{'__rectColor'};
}

sub _rectColor1{
	my $self=shift;
	$self->{'__rectColor1'} = shift // return $self->{'__rectColor1'};
}

sub _rectColor2{
	my $self=shift;
	$self->{'__rectColor2'} = shift // return $self->{'__rectColor2'};
}

sub _font{
	my $self=shift;
	$self->{'__font'} = shift // return $self->{'__font'};
}

sub _fontSize{
	my $self=shift;
	$self->{'__fontSize'} = shift // return $self->{'__fontSize'};
}

sub _fontColor{
	my $self=shift;
	$self->{'__fontColor'} = shift // return $self->{'__fontColor'};
}

sub _labelWidth{
	my $self=shift;
	$self->{'__labelWidth'} = shift // return $self->{'__labelWidth'};
}

sub _mouseText{
	my $self=shift;
	$self->{'__mouseText'} = shift // return $self->{'__mouseText'};
}

sub _printQueryID{
	my $self=shift;
	$self->{'__printQueryID'} = shift // return $self->{'__printQueryID'};
}

sub _mainLabelY{
	my $self=shift;
	$self->{'__mainLabelY'} = shift // return $self->{'__mainLabelY'};
}

sub _minimumPixelWidth{
	my $self=shift;
	$self->{'__minimumPixelWidth'} = shift // return $self->{'__minimumPixelWidth'};
}

sub _spaceFromRightEdge{
	my $self=shift;
	$self->{'__spaceFromRightEdge'} = shift // return $self->{'__spaceFromRightEdge'};
}

sub _scaleY{
	my $self=shift;
	$self->{'__scaleY'} = shift // return $self->{'__scaleY'};
}

sub _scaleLabelY{
	my $self=shift;
	$self->{'__scaleLabelY'} = shift // return $self->{'__scaleLabelY'};
}

sub updateSequences{
	my($self,$qID, $qStart, $qEnd, $qSeqLength, $rID, $rStart, $rEnd, $rSeqLength, $prevQID, $prevRID)=@_;
	
	#get correct sizing
	if($qStart > $qEnd){
		my $tem=$qStart;
		$qStart=$qEnd;
		$qEnd=$tem;
	}
	
	if($rStart > $rEnd){
		my $temp = $rStart;
		$rStart=$rEnd;
		$rEnd=$temp;
	}
	
	if($prevRID ne $rID){
		$self->addToStart($rStart,$qStart);
	}
	
	$self->initializeSequence($qID, $qSeqLength);	
	
}

sub addToStart{
	my($self,$refStart, $queryStart)=@_;
	
	if($queryStart > $refStart){
		
	}
	
}

sub initializeSequence{
	my($self,$query, $seqLen)=@_;
	
	if(!defined $self->_allSequences($query)){
		$self->_allSequences($query,$self->_printCharacter() x $seqLen);
		$self->_maxPosition(0);
	}
	
	
}

sub getSequenceSizesFromFile{
	my($self, $fileToOpen)=@_;
	
	open(FILEIN, "<", $fileToOpen) or die "$!";
	
	local $/=">";
	while(my $aLine=<FILEIN>){
		next if $aLine eq '>';

		my $header;
		my $sequence;
		if($aLine =~ m/^(.*)\n/){
			$header=$1;
			$header =~ s/\W/_/g;
			$aLine =~ s/^(.*)\n//;
			$aLine =~ s/[\W>]//g;
			$sequence=$aLine;
			my $sequenceSize=length($sequence);
			$self->_sequenceSizes->{$header}=$sequenceSize;		
		}
		else{
			print STDERR "could not find fasta file!\n";
			exit(1);
		}		
	}	
	close FILEIN;
}

sub convertSpaces{
	my($self,$stringToConvert)=@_;
	
	$stringToConvert =~ s/\W/_/g;
	return $stringToConvert;
}

sub doneSequence{
	my($self,$doneSeq)=@_;
	
	$self->_sequenceSizes->{$doneSeq}='done';
}

sub getBiggestSequence{
	my($self)=@_;

	my $biggestSize=0;
	my $biggestSizeName='init';
	
	foreach my $seqSizes(keys %{$self->_sequenceSizes}){	
	
		if($self->_sequenceSizes->{$seqSizes} ne 'done'){
			if($self->_sequenceSizes->{$seqSizes} > $biggestSize){
				$biggestSize = $self->_sequenceSizes->{$seqSizes};
				$biggestSizeName=$seqSizes;
			}
		}		
	}	
	
	if($biggestSizeName eq 'init'){
		$self->_biggestSequence('done');
	}
	else{
			$self->_biggestSequence($biggestSizeName);
	}
}

sub getPanGenomeSize{
	my($self)=@_;
	
	my $panSize=0;
	
	foreach my $seqSizes(reverse sort {$a<=>$b} values %{$self->_sequenceSizes()}){	
		$panSize=$panSize + $seqSizes;		
	}	
	$self->_panGenomeSize($panSize);
}




sub getRefSeq{
	my($self,$name)=@_;
	
	my $refValue= FileInteraction::Fasta::SequenceName->new($name);
	return $refValue->name();
}


sub updateValueY{
	my($self, $increment)=@_;

	my $prevRefSeq=$self->getRefSeq($self->_previousQueryID);
	my $currRefSeq=$self->getRefSeq($self->_queryID);
	
	if($prevRefSeq ne $currRefSeq){		
		$self->_newID('true');
		return $increment;
	}
	else{
		$self->_newID('false');
		return 0;
	}	
}

sub passesFilters{
	my($self)=@_;
	
	my $rValue=1;
	
	if($self->_queryLength() < $self->_minimumMatchLength()){
		$rValue=0;
	}	
	
	if($self->_queryID() eq $self->_refID()){
		$rValue=0;
	}
	return $rValue;
}

sub isNucmerLine{
	my($self,$lineToTest)=@_;
	
	if($lineToTest =~ /\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t/){
		return 1;
	}
	else{
		return 0;
	}
}

sub queryLength{
	my($self,$fileName)=@_;
	my $qLength;
	
	open(QFILE,"<", $fileName) or die "$!";
	while(my $line=<QFILE>){
		$line =~ s/[\n\r\f]//g;
		if($.==1){
			my @lineArray=split(/\t/,$line);
			$qLength=$lineArray[8];
		}
		else{
			last;
		}
	}	
	close QFILE;
	$self->_queryLength($qLength);
}

sub getCurrentLine{
	my($self,$line)=@_;
			
	$line =~ s/[\n\r\f\"]//g;
	
	my @lineArray=split(/\t/,$line);
	
	$self->_queryID('init') unless $self->_queryID();
	$self->_refID('init') unless $self->_refID();
	
	my $prevRefID= $self->_refID();
	my $prevQueryID=$self->_queryID();
	
	#ref and query ID get transformed
	#$lineArray[10] =~ s/\W/_/g;
	#$lineArray[9] =~ s/\W/_/g;
	#
	#my $queryName = FileInteraction::Fasta::SequenceName->new($lineArray[10]);
	#my $refName =  FileInteraction::Fasta::SequenceName->new($lineArray[9]);

	my $queryName = $lineArray[10];
	my $refName =  $lineArray[9];
	
	$self->_queryStart($lineArray[2]);
	$self->_queryEnd($lineArray[3]);
	$self->_refStart($lineArray[0]);
	$self->_refEnd($lineArray[1]);
	$self->_refID($refName);
	$self->_previousRefID($prevRefID);
	$self->_previousQueryID($prevQueryID);
	$self->_queryID($queryName);
	$self->_percentID($lineArray[6]);
	$self->_refLength($lineArray[4]);
	$self->_queryLength($lineArray[5]);	
	$self->_refSequenceLength($lineArray[7]);
	$self->_querySequenceLength($lineArray[8]);		
}


sub prettyX{
	my($self)=@_;
	
	if($self->_valueX() < $self->_labelWidth()){
		return $self->_labelWidth();
	}	
	elsif($self->_valueX() > ($self->_graphicWidth() - $self->_spaceFromRightEdge())){
		return ($self->_graphicWidth() - $self->_spaceFromRightEdge());
	}
	else{
		return $self->_valueX();
	}		
	
}

sub createPrintValues{
	my($self,$qID, $qStart, $qEnd, $refStart, $refEnd)=@_;
		
	#make sure qStart < qEnd	
	if($qStart > $qEnd){
		my $tempV=$qStart;
		$qStart = $qEnd;
		$qEnd = $tempV;
	}
	
	#make sure refStart < refEnd	
	if($refStart > $refEnd){
		my $tempV=$refStart;
		$refStart = $refEnd;
		$refEnd = $tempV;
	}
	
	my $length = $refEnd - $refStart + 1;
	my $valuex=int(($refStart / $self->_bpPerUnit())+0.5);
	my $rectWidth=int(($length / $self->_bpPerUnit())+0.5);
		
	#set a minimum pixel width
	if($rectWidth < $self->_minimumPixelWidth()){
		$self->_rectWidth($self->_minimumPixelWidth());
	}
	else{
		$self->_rectWidth($rectWidth);
	}	
	
	$self->_valueX($valuex);

	$self->_printQueryID("$qID" . '(' . "$qStart" . '..' . "$qEnd" . ')');	
}

sub alternateRectColor{
	my($self)=@_;
	
	if($self->_rectColor() eq $self->_rectColor1()){
		$self->_rectColor($self->_rectColor2());
	}
	elsif($self->_rectColor() eq $self->_rectColor2()){
		$self->_rectColor($self->_rectColor1());	
	}
	else{
		print STDERR "color picking problem\n";
		exit(1);
	}	
}

sub getRectWidth{
	my($self,$seqBP)=@_;

	$self->_rectWidth($seqBP / $self->_bpPerUnit());
}

sub bpPerUnit{
	my($self,$sequenceLength,$imageWidth)=@_;
	
	$self->_bpPerUnit($sequenceLength / $imageWidth);
}


1;