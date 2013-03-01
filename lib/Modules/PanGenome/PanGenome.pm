#!/usr/bin/env perl

=pod

=head1 NAME

Modules::PanGenome::PanGenome - 

=head1 SYNOPSIS


	

=head1 DESCRIPTION




=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: 

=head1 AUTHOR

Your name (yourname@email.com)

=head2 Methods

=cut

package Modules::PanGenome::PanGenome;

#includes
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Carp;
use Log::Log4perl;
use Parallel::ForkManager;
use Modules::Alignment::BlastResultFactory;
use Modules::Alignment::SNPFinder;
use Modules::Fasta::MultiFastaSequenceName;
use Role::Tiny::With;

with 'Roles::CombineFilesIntoSingleFile';

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


=head3 _initialize

Initializes the logger.
Assigns all values to class variables.
Anything else that the _initialize function does.

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 

    $self->logger->debug("Logger initialized in Modules::PanGenome::PanGenome");  

    my %params = @_;

    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::PanGenome::PanGenome");
		}
	}	


	#default values
	unless(defined $self->panGenomeOutputFile){
		$self->panGenomeOutputFile($self->outputDirectory . 'pan_genome.txt');
	}

	unless(defined $self->coreSnpsOutputFile){
		$self->coreSnpsOutputFile($self->outputDirectory . 'core_snps.txt');
	}

	$self->_currentResult(0);
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

=head3 xmlFiles

The BLAST-output XML files used for core/accessory processing.

=cut

sub xmlFiles{
	my $self=shift;
	$self->{'_xmlFiles'}=shift // return $self->{'_xmlFiles'};
}

=head3 percentIdentityCutoff

Threshold for sequence identity of a locus to be considered 'core'.
If below, considered 'accessory'

=cut

sub percentIdentityCutoff{
	my $self=shift;
	$self->{'_percentIdentityCutoff'}=shift // return $self->{'_percentIdentityCutoff'};
}

=head3 coreGenomeTheshold

The number of genomes that a locus must be present in to be considered 'core'.
Below this number the locus is considered 'accessory'

=cut

sub coreGenomeThreshold{
	my $self=shift;
	$self->{'_coreGenomeThreshold'}=shift // return $self->{'_coreGenomeThreshold'};
}

=head3 numberOfCores

The number of forks that will be made by Panseq

=cut

sub numberOfCores{
	my $self=shift;
	$self->{'_numberOfCores'}=shift // return $self->{'_numberOfCores'};
}


=head3 resultArraySize

The number of BLAST results that are processed at a time per fork

=cut

sub resultArraySize{
	my $self=shift;
	$self->{'_resultArraySize'}=shift // return $self->{'_resultArraySize'};
}

=head3 outputDirectory

Directory that all files are output to.

=cut

sub outputDirectory{
	my $self=shift;
	$self->{'_outputDirectory'}=shift // return $self->{'_outputDirectory'};
}

=head3 muscleExecutable

Absolute path to the Muscle (http://www.drive5.com/muscle/) alignment program.

=cut

sub muscleExecutable{
	my $self=shift;
	$self->{'_muscleExecutable'}=shift // return $self->{'_muscleExecutable'};
}

=head3 queryFile

All the genomes for the pan-genome are present in this file (Modules::Setup::PanseqFiles)
It is used for creating an ordered list of all names.
This ensures that 'core' / 'accessory' can be properly determined.

=cut

sub queryFile{
	my $self=shift;
	$self->{'_singleQueryFile'}=shift // return $self->{'_singleQueryFile'};
}

=head3 _mfsn 

"multi-fasta sequence name" object.
Stores all of the Modules::Fasta::SequenceName->names for the queryFile

=cut

sub _mfsn{
	my $self=shift;
	$self->{'__mfsn'}=shift // return $self->{'__mfsn'};
}

=head3 accessoryType

Determines whether a binary, percent ID, or actual sequence will be output for the accessory table.

=cut

sub accessoryType{
	my $self=shift;
	$self->{'_accessoryType'}=shift // return $self->{'_accessoryType'};
}


=head3 _orderedNames

An ordered list of names from the _mfsn hash as an array reference.
Computed as a class variable to prevent needless re-computation
every time an array needs to iterate over all names.

=cut

sub _orderedNames{
	my $self=shift;
	$self->{'__orderedNames'}=shift // return $self->{'__orderedNames'};
}


=head3 panGenomeOutputFile

The file that holds the tab-delimited presence / absence, %ID or sequence data
for each genome. This defaults to 'pan_genome.txt' unless overridden.

=cut

sub panGenomeOutputFile{
	my $self=shift;
	$self->{'_panGenomeOutputFile'}=shift // return $self->{'_panGenomeOutputFile'};
}

=head3 coreSnpsOutputFile

The file that holds the tab-delimited SNP data for all "core" regions.
Defaults to 'core_snps.txt' unless overridden.

=cut

sub coreSnpsOutputFile{
	my $self=shift;
	$self->{'_coreSnpsOutputFile'}=shift // return $self->{'_coreSnpsOutputFile'};
}

=head3 _currentResult

Stores the blastXML number across all files and forks.
This number is multiplied by a billion in _processQueue.

=cut

sub _currentResult{
	my $self=shift;
	$self->{'__currentResult'}=shift // return $self->{'__currentResult'};
}


sub run{
	my $self=shift;

	$self->_generateOrderedNamesArray();

	#process all XML files
	foreach my $xml(@{$self->xmlFiles}){
		$self->_processBlastXML($xml);
		unlink $xml;
	}
	$self->_combineTempFiles();
}

=head3 _combineTempFiles

The coreTemp_ and accessoryTemp_ files generated need to be combined into single files.
This does that.
Default filenames are core_snps.txt, pan_genome.txt.
Utilizes Roles::CombineFilesIntoSingleFile
If two or more cores are available, will do both gathers at the same time.

=cut


sub _combineTempFiles{
	my $self = shift;

	my $fileNamesRef = $self->_getFileNamesFromDirectory($self->outputDirectory);
	my $accessoryTempFilesRef = $self->_getAllTempFiles('accessoryTemp_',$fileNamesRef);
	my $coreTempFilesRef = $self->_getAllTempFiles('coreTemp_',$fileNamesRef);

	my $forker = Parallel::ForkManager->new($self->numberOfCores);
	for my $count(1..2){
		$forker->start and next;
			if($count==1){
				$self->_createHeaderLineForCombinedOutputFile($self->panGenomeOutputFile);
				$self->_combineFilesIntoSingleFile(
					$accessoryTempFilesRef,
					$self->panGenomeOutputFile,
					1
				);
				$self->_removeTempFiles($accessoryTempFilesRef);
			}
			elsif($count==2){
				$self->_createHeaderLineForCombinedOutputFile($self->coreSnpsOutputFile);
				$self->_combineFilesIntoSingleFile(
					$coreTempFilesRef,
					$self->coreSnpsOutputFile,
					1
				);
				$self->_removeTempFiles($coreTempFilesRef);
			}else{
				$self->logger->logconfess("Count is $count, but should not be greater than 2");
			}

		$forker->finish;
	}
	$forker->wait_all_children;
}

=head3 _removeTempFiles

After combining the temp files, remove them.
Takes in an arrayRef of file names.

=cut

sub _removeTempFiles{
	my $self=shift;
	my $fileNamesRef=shift;

	foreach my $file(@{$fileNamesRef}){
		#unlink $file;
	}
}


=head3 _createHeaderLineForCombinedOutputFile

The header line consists of tab-delimited genome names over the representative data columns.
Position columns and contigName columns do not have the names redundantly printed, but they are in the same order.

=cut

sub _createHeaderLineForCombinedOutputFile{
	my $self=shift;
	my $outputFile=shift;

	my $outFH = IO::File->new('>' . $outputFile) or die "Could not open $outputFile $!";
	foreach my $name(@{$self->_orderedNames}){
		$outFH->print("\t" . $name);
	}
	$outFH->print("\n");
	$outFH->close();
}


=head3 _getAllTempFiles

Return a list of names that match the first argument passed to the sub.
This allows the gathering of all temp files of a particular type.
The match is conducted against a list of file names passed as the second argument.

=cut

sub _getAllTempFiles{
	my $self=shift;
	my $term = shift;
	my $filesRef=shift;

	my @matchedFiles;
	foreach my $file(@{$filesRef}){
		if($file =~ m/\Q$term\E/){
			push @matchedFiles, $file;
		}
	}
	return \@matchedFiles;
}

sub _generateOrderedNamesArray{
	my $self=shift;

	#generate a hash of all names in the query file
	my $multiFastaSN = Modules::Fasta::MultiFastaSequenceName->new(
		'fileName'=>$self->queryFile
	);
	$self->_mfsn($multiFastaSN);

	#order the hash for a consistent order of names
	$self->_orderedNames([sort keys %{$self->_mfsn->sequenceNameHash}]);
}


sub _processBlastXML {
	my $self = shift;
	my $blastXMLFile    = shift;

	$self->logger->info("Processing Blast XML file $blastXMLFile with " . $self->numberOfCores . ' cores and resultArraySize of ' . $self->resultArraySize);

	my $forker        = Parallel::ForkManager->new($self->numberOfCores);
	my @resultArray;

	#create filehandle
	my $xmlFH = IO::File->new( '<' . $blastXMLFile ) or die "$!";
	my $xmler = Modules::Alignment::BlastResultFactory->new($xmlFH);

	while ( my $result = $xmler->nextResult) {
		push @resultArray, $result;

		if ( scalar(@resultArray) == $self->resultArraySize ) {
			$self->logger->debug("Starting fork. Result " . $self->_currentResult);
			$forker->start and next;
				my ($accRef,$coreRef,$resultNumber)=$self->_processQueue(\@resultArray);
				$self->_emptyBuffers($accRef,$coreRef,$resultNumber);
			$forker->finish;
		}
	}continue {
		$self->_currentResult($self->_currentResult +1);
		if ( scalar(@resultArray) == $self->resultArraySize ) {
			$self->logger->debug("Result array emptied");
			@resultArray = ();
		}
	}
	#$account for stored but unprocessed items
	if(scalar(@resultArray) > 0){
		my ($accRef,$coreRef,$resultNumber)=$self->_processQueue(\@resultArray);
		$self->_emptyBuffers($accRef,$coreRef,$resultNumber);
	}
	
	$forker->wait_all_children();
	$xmlFH->close();
}


=head3 _processQueue

Takes in the BLAST XML data and extracts the core / accessory alignment information from each.
Core is only output if it meets the criteria.
Accessory +/- is output for all loci.
As the number of XML results in the array can vary, the output number is multiplied to give each fork
a large range to work with.

=cut

sub _processQueue {
	my ($self) = shift;

	my $arrayRef        = shift;
	my $locusNumber    = $self->_currentResult;

	my @coreOutputBuffer;
	my @accessoryOutputBuffer;

	#give each temp process a billion number range to work with
	my $resultNumber = $locusNumber * 1000000000;

	foreach my $item ( @{$arrayRef} ) {
		$resultNumber++;
		$self->logger->debug("Processing item $resultNumber"); 

		my $type = $self->_getCoreAccessoryType($item);
		my ($accessoryLine,$coreLines)=$self->_processResult($item, $type, $resultNumber);

		if(defined $accessoryLine){
			push @accessoryOutputBuffer, $accessoryLine;
		}
		else{
			$self->logger->logconfess("accessoryLine should be defined for every input");
		}

		if(defined $coreLines){
			push @coreOutputBuffer, $coreLines;
		}
		else{
			$self->logger->debug("No core lines for $item $resultNumber");
		}
	}
	return(\@accessoryOutputBuffer,\@coreOutputBuffer,$resultNumber);
}


=head3 _emptyBuffers

Takes in an arrayRef to the accessory lines (an array of lines for printing)
and an arrayRef of arrayRefs for core lines (an array where every item is itself an arrayRef)
These need to get printed to temp files in the $self->outputDirectory that will be combined at the end.
The result number ensures no duplicate temp file names.

=cut

sub _emptyBuffers{
	my $self=shift;

	my $accessoryRef = shift;
	my $coreRef = shift;
	my $resultNumber = shift;

	#accessory items have previously been checked (in _processQueue) for definedness
	my $accOutFH = IO::File->new('>>' . $self->outputDirectory . 'accessoryTemp_' . $resultNumber) or $self->logger->logdie("$!");
	foreach my $accLine(@{$accessoryRef}){
		$accOutFH->print($accLine . "\n");
	}
	$accOutFH->close();

	if(defined $coreRef){
		my $coreOutFH = IO::File->new('>>' . $self->outputDirectory . 'coreTemp_' . $resultNumber) or $self->logger->logdie("$!");
		foreach my $arrayRef(@{$coreRef}){
			foreach my $line(@{$arrayRef}){
				unless(defined $line){
					next;
				}
				$coreOutFH->print($line . "\n");
			}
		}
		$coreOutFH->close();
	}
	else{
		$self->logger->logwarn("coreRef buffer undefined in Modules::Alignment::PanGenome::_emptyBuffers resultNumber: $resultNumber");
	}
}

sub _getCoreAccessoryType {
	my ($self) = shift;

	my $result = shift;
	my $returnType;
	my $numberOverSequenceCutoff = 0;

	#if there is no blast result hit, deal with it
	if ( !defined $result->hitHash ) {
		$self->logger->debug("No blast result hit for result");
		$returnType = 'accessory';
	}
	else {
		foreach my $hit ( keys %{ $result->hitHash } ) {
			my $hitObj = $result->hitHash->{$hit};
			my $sequenceCoverage = $self->_getSequenceCoverage( $hitObj, $result->query_len );
			$numberOverSequenceCutoff++ if ( $sequenceCoverage >= $self->percentIdentityCutoff );
		}

		if ( $numberOverSequenceCutoff >= $self->coreGenomeThreshold ) {
			$returnType = 'core';
		}
		else {
			$returnType = 'accessory';
		}
	}
	$self->logger->debug("$returnType number over cutoff: $numberOverSequenceCutoff");
	return $returnType;
}


sub _getSequenceCoverage {
	my $self = shift;

	my $hit         = shift // $self->logger->logconfess('hit is required in Modules::PanGenome::_getSequenceCoverage');
	my $queryLength = shift // $self->logger->logconfess('queryLength is required in Modules::PanGenome::_getSequenceCoverage');
	return ( ( $hit->hsp_align_len - ( $hit->hsp_align_len - $hit->hsp_identity ) ) / $queryLength * 100 );
}

=head3 _processResult

Takes a single BLAST result and extracts both the core and accessory information.
Every result should have an accessory entry (even if all negative)
while a core result is only calculated if the restrictions are met.
In preparation for calculating a core result, an array of fasta sequences is produced,
for direct use in the Muscle alignment program.

=cut

sub _processResult {
	my $self = shift;

	my $result = shift;
	my $type = shift;
	my $resultNumber=shift;

	my $resultHash={};
	my @fastaArray;
	my %startBpHash;
	my $coreResultArrayRef;

	if ( defined $result->hitHash ) {
		foreach my $hit ( keys %{ $result->hitHash } ) {
			my $hitObj = $result->hitHash->{$hit};
			$hitObj->setSequence();    #this ensures start <  end bp

			my $sequenceCoverage = $self->_getSequenceCoverage( $hitObj, $result->query_len );

			$resultHash->{$hit}->{'accessoryValue'} = $self->_getAccessoryValue($hitObj);
			$resultHash->{$hit}->{'fasta'}= $hitObj->hit_def;
			$resultHash->{$hit}->{'startBp'}=$hitObj->hsp_hit_from;

			#create a fasta-file format array if core
			if($type eq 'core'){				
 				push @fastaArray, ( '>' . $hitObj->hit_def . "\n" . $hitObj->hsp_hseq . "\n");
 				$startBpHash{ $hitObj->hit_def } = $hitObj->hsp_hit_from;
			}
		}    #end of foreach
		if($type eq 'core'){
			$coreResultArrayRef = $self->_getCoreResult(\@fastaArray, \%startBpHash, $resultNumber);
		}
	}#end of if
	else{
		$self->logger->info("Result :" . $result->name . " has no hits!");
	}
	my $accessoryResultString = $self->_getAccessoryResult($resultNumber,$resultHash);

	#return the accessory and core results
	return($accessoryResultString,$coreResultArrayRef);
}


=head3 _getAccessoryValue

Returns either a 0, 1, %ID or the sequence for a particular locus in a given genome.
Based on the sequenceCoverage and the accessoryType selected by the user.

=cut

sub _getAccessoryValue{
	my $self=shift;
	my $sequenceCoverage=shift;
	my $hitObj=shift;

	my $value;
	if ( $self->accessoryType eq 'binary' ) {
		if ( $sequenceCoverage >= $self->percentIdentityCutoff ) {
			$value = 1;
		}
		else {
			$value = 0;
		}
	}
	elsif ( $self->accessoryType eq 'percent' ) {
		$value = $sequenceCoverage;
	}
	elsif ( $self->accessoryType eq 'sequence' ) {
		$value = $hitObj->hsp_hseq;
	}
	else {
		$self->logger->fatal("incorrect accessoryType specified!");
	}
	return $value;
}

sub _getAccessoryResult{
	my $self=shift;
	my $resultNumber = shift;
	my $resultHashRef =shift;

	#create the output in correct order
	#the first element of the array is the name of the locus
	#the following are the results per individual genome
	my $counter=0;

	#add query name as first element of array
	my $returnLine='locus_' . $resultNumber;
	my $valueLine='';
	my $positionLine='';
	my $contigNameLine='';

	foreach my $query (@{$self->_orderedNames}) {
		$counter++;
		$valueLine .= "\t";
		$positionLine .="\t";
		$contigNameLine .="\t";

		if ( defined $resultHashRef->{$query} ) {
			$valueLine .= $resultHashRef->{$query}->{'accessoryValue'};
			$positionLine .=$resultHashRef->{$query}->{'startBp'};
			$contigNameLine .=$resultHashRef->{$query}->{'fasta'};
		}
		else {
			$valueLine .= 0;
			$positionLine .='-';
			$contigNameLine .= '-';
		}
	}
	$returnLine .= ($valueLine . $positionLine . $contigNameLine );
	return $returnLine;
}

sub _getCoreResult {
	my $self = shift;

	my $fastaArrayRef = shift;
	my $startBpHashRef = shift;
	my $resultNumber=shift;

 	#create temp files for muscle
 	my $tempInFile = $self->outputDirectory . 'muscleTemp_in' . $resultNumber;
 	my $tempInFH  = IO::File->new('>'. $tempInFile) or die "$!";
 	my $tempOutFile = $self->outputDirectory . 'muscleTemp_out' . $resultNumber;
 	my $tempOutFH = IO::File->new('>' . $tempOutFile) or die "$!";
	
 	$tempInFH->print(@{$fastaArrayRef});	
 	$self->_runMuscle($tempInFile, $tempOutFile);
	
# 	#close the open FH
	$tempInFH->close();
	$tempOutFH->close();
	
# 	#open the output for reading
	my $resultInFH = IO::File->new('<'. $tempOutFile) or die "$!";	
	my @alignedFastaSeqs = $resultInFH->getlines();
	$resultInFH->close();

# 	#delete temp files
	unlink $tempInFile;
	unlink $tempOutFile;

# 	#add SNP information to the return
 	my $snpDetective = Modules::Alignment::SNPFinder->new(
 		'orderedNames'=>$self->_orderedNames,
 		'alignedFastaSequences'=>\@alignedFastaSeqs,
 		'resultNumber'=>$resultNumber,
 		'startBpHashRef'=>$startBpHashRef,
 		'resultNumber'=>$resultNumber
 	);	
 	my $snpDataArrayRef = $snpDetective->findSNPs();
 	#this returns undef if there are no SNPs
 	return $snpDataArrayRef;
}

sub _runMuscle{
	my $self=shift;
	my $inFile = shift;
	my $outFile = shift;
	
	my $systemLine = $self->muscleExecutable . ' -in ' . $inFile . ' -out ' . $outFile . ' -maxiters 3 -quiet';
	system($systemLine);
}

1;





