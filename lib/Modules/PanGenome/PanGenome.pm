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
use diagnostics;
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
	foreach my $xml(sort @{$self->xmlFiles}){
		$counter++;
		$forker->start and next;
			
			$self->_processBlastXML($xml,$counter);
			unlink $xml;
			
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
		$self->settings->baseDirectory . 'coreGenomeFragments.fasta'
	);

	$self->_combineFilesOfType(
		'accessoryGenomeFragments',
		$self->settings->baseDirectory . 'accessoryGenomeFragments.fasta'
	);

	$self->_combineFilesOfType(
		'locus_alleles',
		$self->settings->baseDirectory . 'locus_alleles.fasta'
	);


	#combine separate genome snp/binary phylip strings
	for my $i(1 .. scalar(@{$self->settings->orderedGenomeNames})){
		
		my $snpString = 'snp.phylip_' . $i . '_';
		$self->_combineFilesOfType(
			$snpString,
			$self->settings->baseDirectory . $snpString,
			0 #do not discard first line
		);

		my $binaryString = 'binary.phylip_' . $i . '_';
		$self->_combineFilesOfType(
			$binaryString,
			$self->settings->baseDirectory . $binaryString,
			0 # do not discard first line
		);
	}


	#generate final phylip files
	$self->_combinePhylipFiles('snp.phylip');
	$self->_combinePhylipFiles('binary.phylip');

	#add entries for query segments that have no Blast hits
	if($self->settings->addMissingQuery){
		$self->logger->debug("queryFile specified as " . $self->settings->queryFile);
		$self->_addQueryWithNoBlastHits();
	}
	
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

	$self->logger->info("Processing Blast output file $blastFile, counter: $counter");
	#this should guarantee a unique number for each result of every Panseq run on the same machine
	#allows up to 1000 SNPs per result
	$counter = $self->_getUniqueResultId($counter);
	
	my $blastResult = Modules::Alignment::BlastResults->new($blastFile,$self->settings);
	$self->logger->debug("About to get the first result");
	my $totalResults=0;
	my $totalSeqLength=0;

	my @finalResults;
	while(my $result = $blastResult->getNextResult){
		$totalResults++;		
		

		my @resultKeys = keys %{$result};
		my $numberOfResults = scalar @resultKeys;
		$self->logger->debug("NOR: $numberOfResults");
		unless($numberOfResults > 0){
			$self->logger->warn("No results for query sequence $totalResults, skipping");
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
		
		$totalSeqLength += (length ($result->{$resultKeys[0]}->[0]->[11]));
		
		#this contains any gaps due to the alignment
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

	}#while result

	$self->_printResults($blastFile, \@finalResults);
	$self->logger->info("Total results: $totalResults");
	$self->logger->info("Total base pairs: $totalSeqLength");
}


sub _printResults{
	my $self = shift;
	my $blastFile = shift;
	my $finalResults = shift;

	$blastFile =~ s/\W//g;

	my $binaryFileName = $self->settings->baseDirectory . 'binary_table.txt' . $blastFile;
	my $binaryTableFH = IO::File->new('>>' . $binaryFileName) or die "$!";
	my $snpTableFH = IO::File->new('>>' . $self->settings->baseDirectory . 'snp_table.txt' . $blastFile) or die "$!";
	my $panGenomeFH = IO::File->new('>>' . $self->settings->baseDirectory . 'pan_genome.txt' . $blastFile) or die "$!";
	my $coreSnpsFH = IO::File->new('>>' . $self->settings->baseDirectory . 'core_snps.txt' . $blastFile) or die "$!";
	my $accessoryFH = IO::File->new('>>' . $self->settings->baseDirectory . 'accessoryGenomeFragments.fasta' . $blastFile) or die "$!";
	my $coreFH = IO::File->new('>>' . $self->settings->baseDirectory . 'coreGenomeFragments.fasta' . $blastFile) or die "$!";
	my $allelesFH = IO::File->new('>>' . $self->settings->baseDirectory . 'locus_alleles.fasta' . $blastFile) or die "$!";	

	my @fileHandles = ($binaryTableFH, $snpTableFH, $panGenomeFH, $coreSnpsFH, $accessoryFH, $coreFH, $allelesFH);

	if(-s $binaryFileName < 1){
		#for table files			
		foreach my $genome(@{$self->settings->orderedGenomeNames}){
			$binaryTableFH->print("\t", $genome);
			$snpTableFH->print("\t", $genome);
		};	

		#for pan_genome and core_snps
		my $headerLine = "\tLocusID\tLocusName\tGenome\tAllele\tStart Bp\tEnd Bp\tContig";
		$panGenomeFH->print($headerLine);
		$coreSnpsFH->print($headerLine);
	}
	
	foreach my $finalResult(@{$finalResults}){
		my $locusInformation = $finalResult->{locusInformation};
		my $genomeResults = $finalResult->{genomeResults};

		$self->logger->debug("foreach finalresult: " . Dumper($locusInformation));
		if($self->settings->nameOrId eq 'name'){
			$self->logger->debug("printing: li->name " . Dumper($locusInformation->{name}));
			$binaryTableFH->print("\n", $locusInformation->{name});
		}
		else{
			$binaryTableFH->print("\n", $locusInformation->{id});
		}

		#fragment files
		if($locusInformation->{pan} eq 'core'){
			$coreFH->print('>lcl|', $locusInformation->{id}, '|', $locusInformation->{name}, "\n", $locusInformation->{sequence}, "\n");
		}
		elsif($locusInformation->{pan} eq 'accessory'){
			$accessoryFH->print('>lcl|', $locusInformation->{id}, '|', $locusInformation->{name}, "\n", $locusInformation->{sequence}, "\n");
		}
		else{
			$self->logger->fatal("Unknown pan-genome fragment type!");
			exit(1);
		}

		#alleles file if required
		if($self->settings->storeAlleles){
			$allelesFH->print("\n", 'Locus ', $locusInformation->{name}, "\n");
		}
				
		my $genomeCounter=1;
		my $snpArray=[];
		my @genomesMissingSnps;
		my $numberOfSnps;

		foreach my $genome(@{$self->settings->orderedGenomeNames}){			
			#Table files
			$self->logger->debug("genomes in final result" . @{$self->settings->orderedGenomeNames});
			if(defined $genomeResults->{$genome}->{binary}){
				#alleles file, print out if multiple copies
				if($self->settings->storeAlleles){
					for my $i(0 .. scalar(@{$genomeResults->{$genome}->{binary}}) -1){
						if(defined $genomeResults->{$genome}->{binary}->[$i]->{sequence}){
							my $contigId = $genomeResults->{$genome}->{binary}->[$i]->{contig_id};
							if($i > 1){
								$contigId .= 'a_' . $i;
							}
							$allelesFH->print('>', $contigId, "\n", $genomeResults->{$genome}->{binary}->[$i]->{sequence}, "\n");
						}						
					}					
				}

				$binaryTableFH->print("\t", $genomeResults->{$genome}->{binary}->[0]->{value});

				my $binaryPhylipFH = IO::File->new('>>' . $self->settings->baseDirectory . 'binary.phylip' . '_' . $genomeCounter . '_' . $blastFile ) or die "$!";
				$binaryPhylipFH->print($genomeResults->{$genome}->{binary}->[0]->{value});
				$binaryPhylipFH->close();	

				#pan-genome output			
				$panGenomeFH->print(
					"\n",
					$locusInformation->{id},
					"\t",
					$locusInformation->{name},
					"\t",
					$genome,
					"\t",
					$genomeResults->{$genome}->{binary}->[0]->{value},
					"\t",
					$genomeResults->{$genome}->{binary}->[0]->{start_bp},
					"\t",
					$genomeResults->{$genome}->{binary}->[0]->{end_bp},
					"\t",
					$genomeResults->{$genome}->{binary}->[0]->{contig_id}
				);
			}
			else{
				$self->logger->fatal("No binary result for $genome");
				exit(1);
			}		
			
			if($locusInformation->{pan} eq 'core'){
			
				if(defined $genomeResults->{$genome}->{snp}){
					my $snpLocusCounter=0;
					
					foreach my $snpId(sort keys %{$genomeResults->{$genome}->{snp}}){
						if(!defined $snpArray->[$snpLocusCounter]->[0]){
							$snpArray->[$snpLocusCounter]->[0] = $snpId;
						}
						elsif($snpArray->[$snpLocusCounter]->[0] ne $snpId){
							$self->logger->fatal("SNP IDs do not match!");
							exit(1);
						}				
						$snpArray->[$snpLocusCounter]->[$genomeCounter]=$genomeResults->{$genome}->{snp}->{$snpId}->{value};							
						$snpLocusCounter++;
					}						
				} #defined SNP for genome
			} #end if core
			$genomeCounter++;	
		}#end genome
		#print the SNPs
		$self->_printSnpData(
			snpArray => $snpArray,
			name => $locusInformation->{name},
			genomeResults => $genomeResults,
			snpTableFH => $snpTableFH,
			coreSnpsFH => $coreSnpsFH,
			blastFile => $blastFile
		);
	}

	foreach my $fh(@fileHandles){
		$fh->print("\n");
		$fh->close();
	}	
}


sub _printSnpData{
	my $self = shift;
	my %params = @_;
	
	my $snpArray = $params{snpArray};

	my %snpString;
	foreach my $snp(0..scalar(@{$snpArray})-1){	
		my $snpId = $snpArray->[$snp]->[0];	
		my $snpString = "\n" . $snpId;

		foreach my $genome(1..scalar(@{$self->settings->orderedGenomeNames})){
			my $snpChar = $snpArray->[$snp]->[$genome] // '-';				
			$snpString .= "\t" . $snpChar;

			if(defined $snpString{$genome}){
				$snpString{$genome} .= $snpChar;
			}
			else{
				$snpString{$genome} = $snpChar;
			}

			my $genome = $self->settings->orderedGenomeNames->[$genome -1];
			my $startBp;
			if(defined $params{genomeResults}->{$genome}->{snp}->{$snpId}->{start_bp}){
				$startBp = $params{genomeResults}->{$genome}->{snp}->{$snpId}->{start_bp};
			}
			else{
				$startBp = 0;
			}
			#per SNP printing
			
			$params{coreSnpsFH}->print(
				"\n", 
				$snpId,
				"\t",
				$params{name},
				"\t",
				$genome,
				"\t",
				$snpChar,
				"\t",
				$startBp,
				"\t",
				$startBp,
				"\t",
				$params{genomeResults}->{$genome}->{binary}->[0]->{contig_id}						
			);
		} #foreach
		$params{snpTableFH}->print($snpString);
	}

	foreach my $genome(keys %snpString){
		#phylip file print
		my $snpPhylipFH = IO::File->new('>>' . $self->settings->baseDirectory . 'snp.phylip' . '_' . $genome . '_' . $params{blastFile} ) or die "$!";
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
	my @matchedFiles = sort grep(/\Q$fileNamePart\E/, @{$fileNames});
	$self->logger->debug("Matched files: @matchedFiles");

	my $combiner = Modules::Setup::CombineFilesIntoSingleFile->new();
	$combiner->combineFilesIntoSingleFile(
		\@matchedFiles,
		$outputFile,
		0, #append
		$firstLine #firstLine (1 discards, 0 keeps)
	); 

	#delete original files
	foreach my $file(@matchedFiles){
		unlink $file;
	}
}

sub _combinePhylipFiles{
	my $self = shift;
	my $type = shift;

	my $namer = Modules::Setup::CombineFilesIntoSingleFile->new();
	my $fileNames = $namer->getFileNamesFromDirectory($self->settings->baseDirectory);
	my @matchedFiles = sort grep(/\Q$type\E/, @{$fileNames});

	my $outFH = IO::File->new('>' . $self->settings->baseDirectory . $type) or die "$!";

	my $counter = 1;
	foreach my $phylipFile(@matchedFiles){
		my $inFH = IO::File->new('<' . $phylipFile) or die "$!";
		
		#only a single string per file, therefore all the data is in the first line
		my $fileContent = $inFH->getline();

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
		unlink $phylipFile;
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

	my $conversionFH = IO::File->new('>' . $self->settings->baseDirectory . 'phylip_name_conversion.txt') or die "$!";
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
		
	my $systemLine = $self->settings->muscleExecutable . ' -in ' . $tempInFile . ' -out ' . $tempOutFile . ' -maxiters 3 -quiet';
	system($systemLine);

	#close the open FH
	$tempInFH->close();	
	my @alignedFastaSeqs = $tempOutFH->getlines();
	$tempOutFH->close();

	# #delete temp files
	unlink $tempInFile;
	unlink $tempOutFile;
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





