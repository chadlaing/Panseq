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
use Modules::Fasta::SequenceName;
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

=head3 _mfsn 

"multi-fasta sequence name" object.
Stores all of the Modules::Fasta::SequenceName->names for the queryFile

=cut

sub _mfsn{
	my $self=shift;
	$self->{'__mfsn'}=shift // return $self->{'__mfsn'};
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
	$self->logger->info("Gathering ordered genome names");
	my ($mfsn,$orderedNames)=$self->_generateOrderedNamesArray();
	$self->_mfsn($mfsn);
	$self->_orderedNames($orderedNames);

	my $forker = Parallel::ForkManager->new($self->settings->numberOfCores);
	my $counter=0;
	#process all XML files
	$self->logger->info("Processing Blast output files");
	foreach my $xml(sort @{$self->xmlFiles}){
		$counter++;
		$forker->start and next;
			
			$self->_processBlastXML($xml,$counter);
			#unlink $xml;
			
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



sub _generateOrderedNamesArray{
	my $self=shift;

	#generate a hash of all names in the query file
	my $multiFastaSN = Modules::Fasta::MultiFastaSequenceName->new(
		'fileName'=>$self->queryFile
	);
	
	return ($multiFastaSN,[sort keys %{$multiFastaSN->sequenceNameHash}])
	#order the hash for a consistent order of names
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
		
		$totalSeqLength += (length ($result->{$resultKeys[0]}->[0]->[11]));
		
		my %locusInformation = (
			id=>$counter,
			name=>$result->{$resultKeys[0]}->[0]->[1],
			sequence=>$result->{$resultKeys[0]}->[0]->[11],
			pan=>$coreOrAccessory
		);	
		
		my %genomeResults;
		foreach my $name(@{$self->_orderedNames}){	
			
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
		
		foreach my $name(sort keys %{$result}){
			my $hitNum=0;

			foreach my $hit(@{$result->{$name}}){
				$hitNum++;
				my $contigId = $hit->[0];

				if($hitNum == 1){		
					$genomeResults{$name}->{binary} = [
						{
							contig_id=>$contigId,
							start_bp =>$hit->[2],
							end_bp =>$hit->[3],
							value =>1,
							id=>$counter
						}
					]								
				}
				else{
					push @{$genomeResults{$name}->{binary}},{
						{
							contig_id=>$contigId,
							start_bp =>$hit->[2],
							end_bp =>$hit->[3],
							value =>1,
							id=>$counter
						}
					}
				}				
			
				if($self->settings->storeAlleles){
					push @{$genomeResults{$name}->{alleles}}, {$name . '_a' . $hitNum => $hit->[10]};
				}									
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
			$self->logger->debug("Core, adding to DB");
			my $coreResults = $self->_getCoreResult($result,$msaHash,$counter);	

			foreach my $cResult(@{$coreResults}){
				my $sName = Modules::Fasta::SequenceName->new($cResult->{contig});
				$genomeResults{$sName->name}->{snp}=
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

	my @fileHandles = ($binaryTableFH, $snpTableFH, $panGenomeFH, $coreSnpsFH);

	if(-s $binaryFileName < 1){
		#for table files			
		foreach my $genome(@{$self->_orderedNames()}){
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

		if($self->settings->nameOrId eq 'name'){
			$binaryTableFH->print("\n", $locusInformation->{name});
		}
		else{
			$binaryTableFH->print("\n", $locusInformation->{id});
		}
				
		my $genomeCounter=1;
		my $snpArray=[];
		foreach my $genome(@{$self->_orderedNames()}){			
			#Table files

			if(defined $genomeResults->{$genome}->{binary}){
				$binaryTableFH->print("\t", $genomeResults->{$genome}->{binary}->[0]->{value});
				
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
			
					$coreSnpsFH->print(
						"\n", 
						$snpId,
						"\t",
						$locusInformation->{name},
						"\t",
						$genome,
						"\t",
						$genomeResults->{$genome}->{snp}->{$snpId}->{value},
						"\t",
						$genomeResults->{$genome}->{snp}->{$snpId}->{start_bp},
						"\t",
						$genomeResults->{$genome}->{snp}->{$snpId}->{start_bp},
						"\t",
						$genomeResults->{$genome}->{binary}->[0]->{contig_id}						
					);
				}
			}			
			$genomeCounter++;			
		}#end genome
		#print the SNPs
		foreach my $row(0..scalar(@{$snpArray})-1){
			$snpTableFH->print("\n", $snpArray->[$row]->[0]);
			
			foreach my $column(1..scalar(@{$snpArray->[$row]})-1){
				$snpTableFH->print("\t", $snpArray->[$row]->[$column]);
			}
		}
	}

	foreach my $fh(@fileHandles){
		$fh->print("\n");
		$fh->close();
	}	
}


sub _combineFilesOfType{
	my $self = shift;
	my $fileNamePart = shift;
	my $outputFile = shift;

	$self->logger->warn("Combining files $fileNamePart");
	#use Roles::CombineFilesIntoSingleFile
	my $fileNames = $self->_getFileNamesFromDirectory($self->settings->baseDirectory);
	my @matchedFiles = grep(/$fileNamePart/, @{$fileNames});

	$self->_combineFilesIntoSingleFile(
		\@matchedFiles,
		$outputFile,
		0, #append
		1 #firstLine
	); 

	#delete original files
	foreach my $file(@matchedFiles){
		unlink $file;
	}
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
	foreach my $genome(@{$self->_orderedNames()}){
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
name is based on the Modules::Fasta::SequenceName->name and


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
	#add all the missing genomes as '-'
	foreach my $genome(@{$self->_orderedNames}){
		unless(defined $blastResult->{$genome}){
			$results{'NA_' . $genome} = '-' x $alignmentLength;
		}
	}
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





