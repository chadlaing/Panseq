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
use Modules::Alignment::SNPFinder;
use Modules::Alignment::BlastResults;
use Modules::Setup::CombineFilesIntoSingleFile;
use Data::Dumper;
use File::Basename;

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
    foreach my $key(sort keys %params){
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
		$self->panGenomeOutputFile($self->settings->baseDirectory . 'pan_genome.txt');
	}
	$self->_currentResult(0);
}


sub panGenome{
	my $self=shift;
	$self->{'_panGenome'} = shift // return $self->{'_panGenome'};	
}


sub _printFH{
    my $self = shift;
    $self->{'__printFH'} = shift // return $self->{'__printFH'};   
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



=head3 queryFile

All the genomes for the pan-genome are present in this file (Modules::Setup::PanseqFiles)
It is used for creating an ordered list of all names.
This ensures that 'core' / 'accessory' can be properly determined.

=cut

sub queryFile{
	my $self=shift;
	$self->{'_singleQueryFile'}=shift // return $self->{'_singleQueryFile'};
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

sub settings{
	my $self=shift;
	$self->{'_settings'} = shift // return $self->{'_settings'};	
}


sub _currentResult{
	my $self=shift;
	$self->{'__currentResult'}=shift // return $self->{'__currentResult'};
}

sub _locusId{
	my $self=shift;
	$self->{'_locusId'}=shift // return $self->{'_locusId'};
}


sub run{
	my $self=shift;
	
	$self->logger->info("Analyzing the pan-genome");

	my $forker = Parallel::ForkManager->new($self->settings->numberOfCores);
	my $counter=0;
	#process all XML files
	$self->logger->info("Processing Blast output files");

	#if there is only one file, there will not be a number associated
	#with it. Also we don't need to sort it.
	my @sortedFiles;
	if(scalar(@{$self->xmlFiles}) == 1){
		push @sortedFiles, @{$self->xmlFiles};
	}
	else{
		@sortedFiles = map {$_->[0]}
					   sort {$a->[1] <=> $b->[1]}
					   map {[$_, $self->_getFileNameNumber($_)]}
					   	   @{$self->xmlFiles};
	}
		
	foreach my $xml(@sortedFiles){
		$counter++;
		$forker->start and next;			
			$self->_processBlastXML($xml,$counter);
			unless($self->logger->is_debug()){
				unlink $xml;
			}						
		$forker->finish;
	}
	$forker->wait_all_children;

	#combine into single files
	$self->_combineFilesOfType(
		'snp_table',
		$self->settings->baseDirectory . 'snp_table.txt'
	);

	$self->_combineFilesOfType(
		'binary_table',
		$self->settings->baseDirectory . 'binary_table.txt'
	);

	$self->_combineFilesOfType(
		'pan_genome',
		$self->settings->baseDirectory . 'pan_genome.txt'
	);

	$self->_combineFilesOfType(
		'core_snps',
		$self->settings->baseDirectory . 'core_snps.txt'
	);

	$self->_combineFilesOfType(
		'coreGenomeFragments',
		$self->settings->baseDirectory . 'coreGenomeFragments.fasta',
		0 #do not discard first line
	);

	$self->_combineFilesOfType(
		'accessoryGenomeFragments',
		$self->settings->baseDirectory . 'accessoryGenomeFragments.fasta',
		0 #do not discard first line
	);

	$self->_combineFilesOfType(
		'locus_alleles',
		$self->settings->baseDirectory . 'locus_alleles.fasta',
		0 #do not discard first line
	);


	#combine separate genome snp/binary phylip strings
	for my $i(1 .. scalar(@{$self->settings->orderedGenomeNames})){
		
		my $snpString = 'snp_phylip_' . $i . '_';
		$self->_combineFilesOfType(
			$snpString,
			$self->settings->baseDirectory . $i . '_combined_snp_phylip',
			0 #do not discard first line
		);

		my $binaryString =  'binary_phylip_' . $i . '_';
		$self->_combineFilesOfType(
			$binaryString,
			$self->settings->baseDirectory . $i . '_combined_binary_phylip',
			0 # do not discard first line
		);
	}

	#generate final phylip files
	$self->_combinePhylipFiles('combined_snp_phylip');
	$self->_combinePhylipFiles('combined_binary_phylip');
	
	$self->logger->info("Processing blast output files complete");
	$self->logger->info("Pan-genome generation complete");
}

=head2 _getPanGenomeHashRef

Returns a hash of all "pan" fragments for the run, with the fasta header as key and the sequence as value.

=cut

sub _getPanGenomeHashRef{
	my $self=shift;
	
	my $panFH = IO::File->new('<' . $self->panGenome) or die "$! Could not open " . $self->panGenome . "\n";
	
	my %panHash;
	my $line = $panFH->getline();
	my $sequence='';
	my $header;
	while($line){
		my $nextLine = $panFH->getline();
		
		$line =~ s/\R//g;
		if($line =~ m/^>/){
			$header = $line;
			$header =~ s/>//;
		}
		else{
			$sequence .= $line;
		}		
		
		if((length($sequence) > 0) && (!defined $nextLine || $nextLine =~ m/^>/)){
			$panHash{$header}=$sequence;
		}		
		$line = $nextLine;
	}	
	$panFH->close();
	return \%panHash;
}


=head2 _getUniqueResultId

Given the time in s since the epoch as a starting point, return a 1000 character range of SNPs to work with
for an individual blast result.

=cut

sub _getUniqueResultId{
	my $self=shift;
	my $seed = shift;
	
	return ($seed * 1000 * time());
}

sub _processBlastXML {
	my $self = shift;
	my $blastFile = shift;
	my $counter=shift;

	my $blastFileName = $counter;

	$self->logger->debug("Processing xml $blastFileName");
	$self->_printFH($self->_createPrintFileHandles($blastFileName));
	$self->_printFH->{binaryTableFH}->print("TESY\n\test\tststsatsadtast\nsdtsadtsadt");	
	
	$self->logger->debug("Processing Blast output file $blastFile, counter: $counter");
	#this should guarantee a unique number for each result of every Panseq run on the same machine
	#allows up to 1000 SNPs per result
	$counter = $self->_getUniqueResultId($counter);
	
	my $blastResult = Modules::Alignment::BlastResults->new($blastFile,$self->settings);
	$self->logger->debug("About to get the first result");
	my $totalResults=0;

	my @finalResults;
	while(my $result = $blastResult->getNextResult){
		$totalResults++;		

		my @resultKeys = keys %{$result};
		my $numberOfResults = scalar @resultKeys;
		$self->logger->debug("NOR: $numberOfResults");
		unless($numberOfResults > 0){
			$self->logger->warn("No results for query sequence $totalResults, skipping" . Dumper($result));
			next;
		}		
		$self->logger->debug("There are $numberOfResults results");
		$counter +=1000;			
	
		my $coreOrAccessory;
		if($numberOfResults >= $self->settings->coreGenomeThreshold){
			$coreOrAccessory='core';
		}
		else{
			$coreOrAccessory='accessory';
		}
		$self->logger->debug("coreOrAccrssory: $coreOrAccessory");
		
		#we want the "original" sequence for the locus
		my $queryName = $result->{$resultKeys[0]}->[0]->[1];
		

		#if we are using a query file that has different input sequences,
		#checking the qsn name will auto-vivify empty hash entities
		#which creates downstream troubles
		#if we are using an input query, we cannot guarantee the query will hit
		#so grab the first result and use it for the locus information, rather
		#than using the query name to guarantee the whole perfect match as would
		#be the case in a pan-genome assessment
		my $qsn;
		if(defined $self->settings->queryFile){
			$qsn = $resultKeys[0];
			$self->logger->debug("queryFile defined, Using $qsn as qsn, queryName: $queryName");
		}
		else{
			$qsn = $self->settings->getGenomeNameFromContig($queryName);
			$self->logger->debug("no queryFile, using $qsn as qsn, queryName: $queryName");
		}

		my %locusInformation = (
			id=>$counter,
			name=>$result->{$qsn}->[0]->[1],
			sequence=>$result->{$qsn}->[0]->[11],
			pan=>$coreOrAccessory
		);	
		$self->logger->debug("locusInformation: " . Dumper(%locusInformation));
	
		my %genomeResults;
		foreach my $name(@{$self->settings->orderedGenomeNames}){	
			
			$genomeResults{$name}={
				binary => [
					{contig_id=>"NA",					
					 start_bp=>0,
					 end_bp=>0,
					 value=>0,
					 id=>$counter
					}
				],
				alleles=>[]
			}			
		}
		
		foreach my $name(@resultKeys){
			my $hitNum=1;

			foreach my $hit(@{$result->{$name}}){
				if($hitNum > $self->settings->allelesToKeep){
					last;
				}
				my $contigId = $hit->[0];

				#remove gaps from stored alleles
				$hit->[10] =~ tr/[\-]//d;

				my %binaryHash = (
					contig_id=>$contigId,
					start_bp =>$hit->[2],
					end_bp =>$hit->[3],
					value =>1,
					id=>$counter
				);				

				if($hitNum == 1){		
					$genomeResults{$name}->{binary} = [
						\%binaryHash
					];
					if($self->settings->storeAlleles){
						$genomeResults{$name}->{alleles}=[$hit->[10]];
					}							
				}
				else{
					push @{$genomeResults{$name}->{binary}},\%binaryHash;
					if($self->settings->storeAlleles){
						push @{$genomeResults{$name}->{alleles}}, $hit->[10];
					}
				}				
			
				$hitNum++;							
			}#foreach hit	
		}#foreach name
		
		if($coreOrAccessory eq 'core'){			
			#generate a MSA of all strains that contain sequence
			#if the locus is core, we will send this MSA to the SNPFinder
			my $msaHash = $self->_getHashOfFastaAlignment(
				$self->_getMsa($result,$counter),
				$result
			);
			
			#'outfmt'=>'"6 
			# [0]sseqid 
			# [1]qseqid 
			# [2]sstart 
			# [3]send 
			# [4]qstart 
			# [5]qend 
			# [6]slen 
			# [7]qlen 
			# [8]pident 
			# [9]length"',
			# [10]sseq,
			# [11]qseq		
	
			#if it is a core result, send to SNP finding
			my $coreResults = $self->_getCoreResult($result,$msaHash,$counter);	

			foreach my $cResult(@{$coreResults}){
				#$self->logger->warn("cResult contig: " . $cResult->{contig});
				my $sName = $self->settings->getGenomeNameFromContig($cResult->{contig}) // $self->logger->logdie("Could not find genome for " . $cResult->{contig});
				#my $sName = $cResult->{contig};
				$genomeResults{$sName}->{snp}=
					{
						$cResult->{locusId}=>{
							start_bp=>$cResult->{startBp},
							value=>$cResult->{value}
						}
					};
			} #foreach cResult
		}#if core

		my %result=(
			locusInformation=>\%locusInformation,
			genomeResults=>\%genomeResults
		);
		push @finalResults, \%result;

		if(scalar(@finalResults) == $self->settings->maxNumberResultsInMemory){
			$self->_printResults($blastFileName,\@finalResults);
			@finalResults=();
		}

	}#while result
	$self->_printResults($blastFileName, \@finalResults);

	my $conversionFile = $self->settings->baseDirectory . 'phylip_name_conversion.txt';
	unless(-e $conversionFile){
		$self->_printConversionInformation($conversionFile);
	}
	$self->logger->info("Total results: $totalResults");
	$self->_closePrintFileHandles();
}


sub _getNameOrId{
	my $self=shift;
	my $locusInformation = shift;

	if($self->settings->nameOrId eq 'name'){
		return("\n" . $locusInformation->{name});
	}
	else{
		return("\n" . $locusInformation->{id});
	}
}


sub _printPanGenomeFastaFiles{
	my $self = shift;
	my $locusInformation = shift;

	#fragment files
	if($locusInformation->{pan} eq 'core'){
		$self->_printFH->{coreFH}->print('>lcl|' 
						 .$locusInformation->{id} 
						 . '|' 
						 . $locusInformation->{name} 
						 . "\n" 
						 . $locusInformation->{sequence} 
						 . "\n");
	}
	elsif($locusInformation->{pan} eq 'accessory'){
		$self->_printFH->{accessoryFH}->print('>lcl|'
											 . $locusInformation->{id} 
											 . '|'
											 . $locusInformation->{name} 
											 . "\n" 
											 . $locusInformation->{sequence}
											 . "\n");
	}
	else{
		$self->logger->fatal("Unknown pan-genome fragment type! " . $locusInformation->{pan});
		exit(1);
	}
}


sub _printAlleleData{
	my $self = shift;
	my $finalResult = shift;
	my $alleleData = shift;
	my $firstFlag =shift;
	my $contigId = shift;

	$self->logger->debug("Print allele data: contig: $contigId");
	#alleles file if required
	if($firstFlag){
		$self->_printFH->{allelesFH}->print("\n" . 'Locus ' . $finalResult->{locusInformation}->{name} . "\n");
	}	
	
	$self->_printFH->{allelesFH}->print('>' . $contigId . "\n" .  $alleleData . "\n");
	return 0;
}


sub _printBinaryTableData{
	my $self = shift;
	my $nameOrId = shift;
	my $data = shift;
	$self->_printFH->{binaryTableFH}->print("\t" . $data);
}


sub _printBinaryPhylipData{
	my $self = shift;
	my $fileName = shift;
	my $data = shift;

	my $binaryPhylipFH = IO::File->new('>>' . $fileName) or die "$!";
	$binaryPhylipFH->print($data);
	$binaryPhylipFH->close();	
}


sub _addToSnpArray{
	my $self = shift;
	my $snpArray = shift;
	my $singleGenome = shift;
	my $genomeCounter = shift;

	if(defined $singleGenome->{snp}){
		my $snpLocusCounter=0;
		
		foreach my $snpId(sort keys %{$singleGenome->{snp}}){
			if(!defined $snpArray->[$snpLocusCounter]->[0]){
				$snpArray->[$snpLocusCounter]->[0] = $snpId;
			}
			elsif($snpArray->[$snpLocusCounter]->[0] ne $snpId){
				$self->logger->fatal("SNP IDs do not match!");
				exit(1);
			}				
			$snpArray->[$snpLocusCounter]->[$genomeCounter]=$singleGenome->{snp}->{$snpId}->{value};							
			$snpLocusCounter++;
		}						
	}
	return $snpArray;
}


sub _getContigId{
	my $self = shift;
	my $contigId = shift;
	my $i = shift;
	
	if($i > 0){
		$contigId .= '_a' . ($i +1);
	}					
	return $contigId;									
}

sub _printPanGenomeData{
	my $self=shift;
	my $locusInformation = shift;
	my $singleGenome = shift;
	my $genome = shift;
		
	$self->_printFH->{panGenomeFH}->print(
		"\n" .
		$locusInformation->{id} .
		"\t" .
		$locusInformation->{name} .
		"\t" .
		$genome .
		"\t" .
		$singleGenome->{value} .
		"\t" .
		$singleGenome->{start_bp} .
		"\t" .
		$singleGenome->{end_bp} .
		"\t" .
		$singleGenome->{contig_id}
	);		
}


sub _printResults{
	my $self = shift;
	my $blastFile = shift;
	my $finalResults = shift;
	
	foreach my $finalResult(@{$finalResults}){
		$self->_printPanGenomeFastaFiles($finalResult->{locusInformation});		
				
		#my $snpArray=[];
		my $genomeCounter = 1;
		my $firstFlag = 1;
		foreach my $genome(@{$self->settings->orderedGenomeNames}){	
			#do all of the single allele printing first
			$self->_printBinaryTableData($self->_getNameOrId($finalResult->{locusInformation})										
										,$finalResult->{genomeResults}->{$genome}->{binary}->[0]->{value});

			$self->_printBinaryPhylipData($self->settings->baseDirectory . $blastFile . '_binary_phylip_' . $genomeCounter . '_' 
										 ,$finalResult->{genomeResults}->{$genome}->{binary}->[0]->{value});

			# if($finalResult->{locusInformation}->{pan} eq 'core'){
			# 	$snpArray = $self->_addToSnpArray($snpArray
			# 									 ,$finalResult->{genomeResults}->{$genome}
			# 									 ,$genomeCounter);	
			# }

			#now we are concerned about the multiple alleles
			#print with them in mind
			if(defined $finalResult->{genomeResults}->{$genome}->{binary}){
				for my $i(0 .. scalar(@{$finalResult->{genomeResults}->{$genome}->{binary}}) -1){
					my $contigId = $self->_getContigId($finalResult->{genomeResults}->{$genome}->{binary}->[$i]->{contig_id}, $i);

					#if the contig is NA, there is no allele data to print
					if($self->settings->storeAlleles && $contigId ne 'NA'){
						$firstFlag = $self->_printAlleleData($finalResult
															,$finalResult->{genomeResults}->{$genome}->{alleles}->[$i]
															,$firstFlag
															,$contigId);
					}			
				
					$self->_printPanGenomeData($finalResult->{locusInformation}
											  ,$finalResult->{genomeResults}->{$genome}->{binary}->[$i]
											  ,$genome);
				} #end allele	
			}#end genomeResult
			$genomeCounter++;
		} #end genome
		#print the SNPs
		$self->_printSnpData(
			#snpArray => $snpArray,
			name => $finalResult->{locusInformation}->{name},
			genomeResults => $finalResult->{genomeResults},
			blastFile => $blastFile
		);
	}#end finalResults		
}


sub _createPrintFileHandles{
	my $self = shift;
	my $blastFile = shift;

	my %allFH = (
		binaryTableFH => (IO::File->new('>' . $self->settings->baseDirectory . $blastFile . '_binary_table') or die "$!"),
		snpTableFH => (IO::File->new('>' . $self->settings->baseDirectory . $blastFile . '_snp_table') or die "$!"),
		panGenomeFH => (IO::File->new('>' . $self->settings->baseDirectory . $blastFile . '_pan_genome') or die "$!"),
		coreSnpsFH => (IO::File->new('>' . $self->settings->baseDirectory . $blastFile . '_core_snps') or die "$!"),
		accessoryFH => (IO::File->new('>' . $self->settings->baseDirectory . $blastFile . '_accessoryGenomeFragments') or die "$!"),
		coreFH => (IO::File->new('>' . $self->settings->baseDirectory . $blastFile . '_coreGenomeFragments') or die "$!"),
		allelesFH => (IO::File->new('>' . $self->settings->baseDirectory . $blastFile . '_locus_alleles') or die "$!")
	);

	#Initialize the table files with headers			
	foreach my $genome(@{$self->settings->orderedGenomeNames}){
		$allFH{binaryTableFH}->print("\t", $genome);
		$allFH{snpTableFH}->print("\t", $genome);
	};	

	#for pan_genome and core_snps
	my $headerLine = "LocusID\tLocusName\tGenome\tAllele\tStart Bp\tEnd Bp\tContig";
	$allFH{panGenomeFH}->print($headerLine);
	$allFH{coreSnpsFH}->print($headerLine);

	return \%allFH;
}



sub _closePrintFileHandles{
	my $self = shift;

	foreach my $fh(keys %{$self->_printFH}){
		$self->_printFH->{$fh}->close();
	}
}


sub _printSnpData{
	my $self = shift;
	my %params = @_;
	
	my $snpArray = $params{snpArray};

	$self->logger->debug("Printing snp data:");
	$self->logger->debug(Dumper(%params));

	my %snpString;
	for my $i(0..scalar(@{$snpArray})-1){	
		my $snpId = $snpArray->[$i]->[0];	
		my $snpString = "\n" . $snpId;

		for my $j(1..scalar(@{$self->settings->orderedGenomeNames})){
			my $snpChar = $snpArray->[$i]->[$j] // '-';				
			$snpString .= "\t" . $snpChar;

			if(defined $snpString{$j}){
				$snpString{$j} .= $snpChar;
			}
			else{
				$snpString{$j} = $snpChar;
			}

			my $genome = $self->settings->orderedGenomeNames->[$j -1];
			
			my $startBp = 0;
			if(defined $params{genomeResults}->{$genome}->{snp}->{$snpId}->{start_bp}){
				$startBp = $params{genomeResults}->{$genome}->{snp}->{$snpId}->{start_bp};
			}
			#per SNP printing
			
			$self->_printFH->{coreSnpsFH}->print(
				"\n" . 
				$snpId .
				"\t" .
				$params{name} .
				"\t" .
				$genome .
				"\t" .
				$snpChar .
				"\t" .
				$startBp .
				"\t" .
				$startBp .
				"\t" .
				$params{genomeResults}->{$genome}->{binary}->[0]->{contig_id}						
			);
		} #foreach
		$self->_printFH->{snpTableFH}->print($snpString);
	}

	#these keys do not need to be sorted, as they print to the file given by the genome number
	#it doesn't matter what order they are printed in
	foreach my $genome(keys %snpString){
		#phylip file print
		my $snpPhylipFH = IO::File->new('>>' . $self->settings->baseDirectory . $params{blastFile} . '_snp_phylip_' . $genome . '_') or die "$!";
		$snpPhylipFH->print($snpString{$genome});
		$snpPhylipFH->close();
	}
}


sub _combineFilesOfType{
	my $self = shift;
	my $fileNamePart = shift;
	my $outputFile = shift;
	my $firstLine = shift // 1;

	my $namer = Modules::Setup::CombineFilesIntoSingleFile->new();
	my $fileNames = $namer->getFileNamesFromDirectory($self->settings->baseDirectory);
	my @matchedFiles = map {$_->[0]}
					   sort {$a->[1] <=> $b->[1]}
					   map {[$_, $self->_getFileNameNumber($_)]}
					   	   grep(/\Q$fileNamePart\E/, @{$fileNames});
	$self->logger->debug("combineFilesOfType Matched files: @matchedFiles");

	my $combiner = Modules::Setup::CombineFilesIntoSingleFile->new();
	$combiner->combineFilesIntoSingleFile(
		\@matchedFiles,
		$outputFile,
		0, #append
		$firstLine #firstLine (1 discards, 0 keeps)
	); 

	#delete original files
	unless($self->logger->is_debug()){
		foreach my $file(@matchedFiles){
			unlink $file;
		}
	}
}


sub _getFileNameNumber{
	my $self = shift;
	my $fileNameWithDir = shift;

	#with File::Basename
	my $fileName = fileparse($fileNameWithDir);

	if($fileName =~ m/^(\d+)_/){
		return $1;
	}
	else{
		$self->logger->fatal("Could not find filename number in file $fileName");
		exit(1);
	}
}

sub _combinePhylipFiles{
	my $self = shift;
	my $type = shift;

	my $namer = Modules::Setup::CombineFilesIntoSingleFile->new();
	my $fileNames = $namer->getFileNamesFromDirectory($self->settings->baseDirectory);
	my @matchedFiles = map {$_->[0]}
					   sort {$a->[1] <=> $b->[1]}
					   map {[$_, $self->_getFileNameNumber($_)]}
					   	  grep(/\Q$type\E/, @{$fileNames});

	$self->logger->debug("matchedFiles: @matchedFiles");

	my $fileName;
	if($type eq 'combined_binary_phylip'){
		$fileName = 'binary.phylip';
	}
	elsif($type eq 'combined_snp_phylip'){
		$fileName = 'snp.phylip';
	}
	else{
		$self->logger->fatal("Unknown phylip file type");
		exit(1);
	}


	my $outFH = IO::File->new('>' . $self->settings->baseDirectory . $fileName) or die "$!";

	my $counter = 1;
	foreach my $phylipFile(@matchedFiles){
		my $inFH = IO::File->new('<' . $phylipFile) or die "$!";
		
		#only a single string per file, therefore all the data is in the first line
		my @fileContents = $inFH->getlines();
		my $fileContent = join('', @fileContents);

		unless(defined $fileContent){
			if($type eq 'binary'){
				$self->logger->logconfess("Binary file is empty. Fatal error");
			}
			else{
				$self->logger->warn("SNP file empty");
				next;
			}
		}
		$fileContent =~ s/\R//g;

		if($counter == 1){
			#number of genome
			$outFH->print(scalar(@matchedFiles), ' ', length($fileContent), "\n");
		}
		$outFH->print($counter, $self->_numberOfSpacesToAdd($counter), ' ', $fileContent, "\n");
		$counter++;
		$inFH->close();

		unless($self->logger->is_debug()){
			unlink $phylipFile;
		}		
	}
	$outFH->close();
}


=head2 _numberOfSpacesToAdd

Phylip requires the name to be exactly 10 characters, including spaces.
This is that.

=cut

sub _numberOfSpacesToAdd{
	my $self=shift;
	my $counter = shift;
	
	my $numberOfSpaces = 10 - length($counter);
	return(" " x $numberOfSpaces);	
}


=head2 _printConversionInformation

Phylip format is limited to a 10-character name field.
In printing the Phylip format, we substitute numbers for names.
This creates a tab-delimited table that lists the conversion information.

=cut

sub _printConversionInformation{
	my $self=shift;	
	my $conversionFile = shift;

	my $conversionFH = IO::File->new('>' . $conversionFile) or die "Could not open $conversionFile$!";
	$conversionFH->print(
		'Number' . "\t" . 'Name' . "\n"
	);

	my $counter=1;
	foreach my $genome(@{$self->settings->orderedGenomeNames}){
		$conversionFH->print($counter . "\t" . $genome . "\n");
		$counter++;
	}
	$conversionFH->close();	
}


=head2 _getMsa

Take in a Modules::Alignment::BlastResult and generate a MSA of all the TOP hit sequences.

=cut

sub _getMsa{
	my $self=shift;
	my $result=shift;	
	my $resultNumber=shift;

	#create temp files for muscle
	my $tempInFile = $self->settings->baseDirectory . 'muscleTemp_in' . $resultNumber;
	my $tempInFH = IO::File->new('>'. $tempInFile) or die "$!";
	my $tempOutFile = $self->settings->baseDirectory . 'muscleTemp_out' . $resultNumber;
	my $tempOutFH = IO::File->new('+>' . $tempOutFile) or die "$!";
	
	#'outfmt'=>'"6 
	# [0]sseqid 
	# [1]qseqid 
	# [2]sstart 
	# [3]send 
	# [4]qstart 
	# [5]qend 
	# [6]slen 
	# [7]qlen 
	# [8]pident 
	# [9]length"',
	# [10]sseq,
	# [11]qseq
	foreach my $resKey(keys %{$result}){
		$tempInFH->print('>' . $result->{$resKey}->[0]->[0] . "\n" . $result->{$resKey}->[0]->[10] . "\n");
	}
		
	my $systemLine = $self->settings->muscleExecutable . ' -in ' . $tempInFile . ' -out ' . $tempOutFile . ' -maxiters 2 -diags -quiet';
	system($systemLine);

	#close the open FH
	$tempInFH->close();	
	my @alignedFastaSeqs = $tempOutFH->getlines();
	$tempOutFH->close();

	# #delete temp files
	unless($self->logger->is_debug()){
		unlink $tempInFile;
		unlink $tempOutFile;
	}
	return \@alignedFastaSeqs;
}

=head2 _getHashOfFastaAlignment

Given the FASTA alignment produced by Muscle, create a hash where the
name is the contig

=cut

sub _getHashOfFastaAlignment{
	my $self = shift;
	my $alignedFastaSequences=shift;
	my $blastResult = shift;
	
	my %results;
	my $header;	
	
	foreach my $line(@{$alignedFastaSequences}){
		$line =~ s/\R//g;
		
		if($line =~ /^>(.+)/){		
			$line =~ s/>//;			
			$header= $line;			
		}
		else{
			if(defined $results{$header}){
				$results{$header} .= $line;
				
			}
			else{
				$results{$header} = $line;
			}					
		}
	}
	my @gn = keys %results;
	my $alignmentLength = length($results{$gn[0]});
	
	return (\%results);
}


sub _getCoreResult {
	my $self = shift;
	my $result=shift;
	my $msaHash=shift;
	my $resultNumber=shift;
	
	#'outfmt'=>'"6 
	# [0]sseqid 
	# [1]qseqid 
	# [2]sstart 
	# [3]send 
	# [4]qstart 
	# [5]qend 
	# [6]slen 
	# [7]qlen 
	# [8]pident 
	# [9]length"',
	# [10]sseq,
	# [11]qseq
	my %startBpHash;
	foreach my $resKey(keys %{$result}){
		$startBpHash{$result->{$resKey}->[0]->[0]}=$result->{$resKey}->[0]->[4];
	}
	
	#add SNP information to the return
	my $snpDetective = Modules::Alignment::SNPFinder->new(
		 'alignedFastaHash'=>$msaHash,
		 'startBpHashRef'=>\%startBpHash,
		 'resultNumber'=>$resultNumber,
		 'frameshift'=>$self->settings->frameshift
	 );	
	 my $snpDataArrayRef = $snpDetective->findSNPs();
	 #this returns undef if there are no SNPs
	 return $snpDataArrayRef;	
}

1;





