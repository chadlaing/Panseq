#!/usr/bin/perl


#This needs some serious modernization, but it works
#Alas
#It also needs some serious documentation, but it works
#Alas

#_fastaQueryFile should be a "pan-genome" of contigs that the reference bar on the image can be based on
#_inputFile refers to a nucmer out.file processed vis show-coords with -l -q -T options. The comparison should be between the
#pan-genome file used as fastaQueryFile but with all the fasta headers removed (eg. a single concatenated sequence with a single fasta header)
#an arbitrary ">panGenome" works fine. 
#eg. ~/MUMmer3.23/nucmer pangenome_single.fasta queryFilesCombined.fasta
# ~/MUMmer3.23/show-coords out.delta -l -q -T > out.coords


use strict;
use warnings;
use diagnostics;
use FindBin;
use lib "$FindBin::Bin/../";
use Visual::SvgDrawer;
use XML::Writer;


my $filer=Visual::SvgDrawer->new();

$filer->_inputFile($ARGV[0]);
$filer->_fastaQueryFile($ARGV[1]);
$filer->_minimumMatchLength($ARGV[2]);

my $x=XML::Writer->new();
open(NUCMERIN, "<", $filer->_inputFile()) or die "$!";

my $imager=Visual::SvgDrawer->new();
$imager->_imageUnits('px');
$imager->_imageHeight(1000);
$imager->_imageWidth(1400);
$imager->_labelWidth(200);
$imager->_graphicWidth($imager->_imageWidth() - ($imager->_labelWidth()*2));
$imager->_valueY(20);
$imager->_valueX($imager->_labelWidth());
$imager->_incrementY(10);
$imager->_rectHeight(20);
$imager->_rectColor('#d70d0d');
$imager->_rectColor1('#fb1212');
$imager->_rectColor2('#d70d0d');
$imager->_font("serif");
$imager->_fontSize('12');
$imager->_fontColor('black');
$imager->_mainLabelY(15);
$imager->_spaceFromRightEdge('800');
$imager->_minimumPixelWidth(1);
$imager->_scaleY(50);
$imager->_scaleLabelY(70);


$x->startTag('svg', 'id'=>'body','overflow'=>'auto', 'width'=>$imager->_imageWidth() . $imager->_imageUnits(),'height'=>$imager->_imageHeight(). $imager->_imageUnits(), 'version'=>'1.1', 'xmlns'=>"http://www.w3.org/2000/svg", 'xmlns:xlink'=>"http://www.w3.org/1999/xlink");
	#add style-sheet properties
	$x->startTag('style', 'type'=>'text/css');
		$x->characters('			
		      <![CDATA[		
		     		svg:root{
		     			overflow: auto;
		     		}
		      ]]>				
		');
	$x->endTag('style');
	
	$x->startTag('title');
		$x->characters('Sequence Alignment');
	$x->endTag('title');
	
	$x->startTag('desc');
		$x->characters('Panseq sequence alignment');
	$x->endTag('desc');
	
	#go through queryfile fasta file
	my $seqAn=Visual::SvgDrawer->new();
	$seqAn->getSequenceSizesFromFile($filer->_fastaQueryFile()); #in _sequenceSizes
	
	$seqAn->getPanGenomeSize(); #save panGenome size to _panGenomeSize
	$imager->bpPerUnit($seqAn->_panGenomeSize(),$imager->_graphicWidth()); #get bp/unit in _bpPerPixel
	
	$x->startTag('text','font-family'=>$imager->_font(), 'font-size'=>$imager->_fontSize(),'fill'=>$imager->_fontColor(), 'x'=>"1", 'y'=>"40");
		$x->characters('Complete Sequence');		
	$x->endTag('text');
	
	$seqAn->getBiggestSequence();
	while($seqAn->_biggestSequence() ne 'done'){

		$imager->getRectWidth($seqAn->_sequenceSizes->{$seqAn->_biggestSequence}); #store in _rectWidth
		
		$x->startTag('rect','id'=>$seqAn->_biggestSequence(), 'x'=>$imager->_valueX(),'y'=>$imager->_valueY(),'width'=>$imager->_rectWidth(),'height'=>$imager->_rectHeight(),'fill'=> $imager->_rectColor());
			$x->startTag('set','attributeName'=>'fill', 'to'=>'black', 'begin'=>'mouseover', 'end'=>'mouseout');
			$x->endTag('set');			
		$x->endTag('rect');
		
		$x->startTag('text','font-family'=>$imager->_font(), 'visibility'=>'hidden', 'font-size'=>$imager->_fontSize(),'fill'=>$imager->_fontColor(), 'x'=>$imager->prettyX(), 'y'=>$imager->_mainLabelY());
			$x->characters($seqAn->_biggestSequence());	
			$x->startTag('set','attributeName'=>'visibility', 'to'=>'visible', 'begin'=>$seqAn->_biggestSequence(). ".mouseover", 'end'=>$seqAn->_biggestSequence() . ".mouseout");
			$x->endTag('set');	
		$x->endTag('text');
		
		$imager->_valueX($imager->_valueX() + $imager->_rectWidth()); #move x value to new position
		$seqAn->doneSequence($seqAn->_biggestSequence());
		$seqAn->getBiggestSequence();
		$imager->alternateRectColor();
		#$seqAn->_biggestSequence('done'); #for single sequence testing
	}
		
	#add a scale to the image
	$x->startTag('text','font-family'=>$imager->_font(),'font-size'=>$imager->_fontSize(),'fill'=>$imager->_fontColor(),'x'=>'0','y'=>$imager->_scaleLabelY());
		$x->characters('Size (bp)');
	$x->endTag('text');
	
	my $width=($imager->_graphicWidth())/10;
	$x->startTag('defs');
		$x->startTag('polyline', 'id'=>'scalePart','fill'=>'none','stroke'=>'black','stroke-width'=>'1',
			'points'=>"0,0 $width,0 $width,10");
		$x->endTag('polyline');		
	$x->endTag('defs');	
	
	my $total=$imager->_labelWidth();
	for(my $piece=1; $piece <=10; $piece++){			
		$x->startTag('use','xlink:href'=>'#scalePart','transform'=>"translate($total," . $imager->_scaleY() .")");	
		$x->endTag('use');
		
		$x->startTag('text','font-family'=>$imager->_font(),'font-size'=>$imager->_fontSize(),'fill'=>$imager->_fontColor(),'x'=>$total,'y'=>$imager->_scaleLabelY());
			$x->characters(int($imager->_bpPerUnit() * ($total - $imager->_labelWidth())));
		$x->endTag('text');
		$total=$total + $width;
	}
	$x->startTag('text','font-family'=>$imager->_font(),'font-size'=>$imager->_fontSize(),'fill'=>$imager->_fontColor(),'x'=>$total,'y'=>$imager->_scaleLabelY());
			$x->characters(int($imager->_bpPerUnit() * ($total - $imager->_labelWidth())));
	$x->endTag('text');
	#END of scale
		
	my $seqImager = Visual::SvgDrawer->new();
	$seqImager->_rectHeight(20);
	$seqImager->_rectColor('#222222');
	$seqImager->_rectColor1('#222222');
	$seqImager->_rectColor2('#333333');
	$seqImager->_font("serif");
	$seqImager->_fontSize('12');
	$seqImager->_fontColor('black');
	$seqImager->_valueY('80');
	$seqImager->_incrementY('30');
	$seqImager->_bpPerUnit($imager->_bpPerUnit());
	$seqImager->_minimumPixelWidth(1);
	$seqImager->_spaceFromRightEdge($imager->_spaceFromRightEdge());
	$seqImager->_imageWidth($imager->_imageWidth());
	$seqImager->_graphicWidth($imager->_graphicWidth());
	$seqImager->_labelWidth($imager->_labelWidth());
	$seqImager->_valueX(0);
	
	while(my $inline=<NUCMERIN>){
		next unless $filer->isNucmerLine($inline); #make sure properly formatted line for analysis
		$filer->getCurrentLine($inline); #updates all file values 				
		$seqImager->createPrintValues($filer->_queryID(),$filer->_queryStart(),$filer->_queryEnd(), $filer->_refStart(), $filer->_refEnd()); #creates _printQueryID, _rectWidth, _valueX
		$seqImager->_valueY($seqImager->_valueY() + $filer->updateValueY($seqImager->_incrementY())); #increments _valueY if different sequence names, sets '_newID' to true if so
		
		
		#add tag if new sequence
		if($filer->_newID() eq 'true'){
			$x->startTag('text', 'font-family'=>$seqImager->_font(), 'font-size'=>$seqImager->_fontSize(), 'fill'=>$seqImager->_fontColor(),'x'=>'0', 'y'=>$seqImager->_valueY()+$seqImager->_fontSize());
				$x->characters($filer->getRefSeq($filer->_queryID()));
			$x->endTag('text')
		}
		#
		
		#even if match is too small, need to print out the new ID of a new strain
		next unless $filer->passesFilters(); #if the match is too small etc., skip
		
		$x->startTag('rect','id'=>'item'."$.", 'x'=>$seqImager->_valueX() + $imager->_labelWidth(),'y'=>$seqImager->_valueY(),'fill'=>$seqImager->_rectColor(), 'width'=>$seqImager->_rectWidth(),'height'=>$seqImager->_rectHeight());
			$x->startTag('set','attributeName'=>'fill', 'to'=>'blue', 'begin'=>'mouseover', 'end'=>'mouseout');
			$x->endTag('set');	
		$x->endTag('rect');	
		
		$x->startTag('text','font-family'=>$seqImager->_font(), 'visibility'=>'hidden', 'font-size'=>$seqImager->_fontSize(),'fill'=>$seqImager->_fontColor(), 'x'=>$seqImager->prettyX(), 'y'=>$imager->_mainLabelY());
			$x->characters($seqImager->_printQueryID());	
			$x->startTag('set','attributeName'=>'visibility', 'to'=>'visible', 'begin'=>"item$.". ".mouseover", 'end'=>"item$." . ".mouseout");
			$x->endTag('set');	
		$x->endTag('text');
	}
		
$x->endTag('svg');

close NUCMERIN;