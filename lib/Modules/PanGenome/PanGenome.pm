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
use DBI;
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

	#construct sqlite DB
	$self->_sqlString({});
	$self->_sqlSelectNumber({});
	$self->_initDb();

	#default values
	unless(defined $self->panGenomeOutputFile){
		$self->panGenomeOutputFile($self->settings->baseDirectory . 'pan_genome.txt');
	}

	$self->_currentResult(0);
}

sub _initDb{
	my $self = shift;
	
	$self->logger->info("Initializing SQL DB");
	#define SQLite db
	my $dbh = (DBI->connect("dbi:SQLite:dbname=" . $self->settings->baseDirectory . "temp_sql.db","","")) or $self->logdie("Could not connect to SQLite DB");
	
	$dbh->do("DROP TABLE IF EXISTS results");
	$dbh->do("DROP TABLE IF EXISTS strain");
	$dbh->do("DROP TABLE IF EXISTS contig");
	$dbh->do("DROP TABLE IF EXISTS locus");
	$dbh->do("DROP TABLE IF EXISTS allele");

	$dbh->do("CREATE TABLE strain(id INTEGER PRIMARY KEY not NULL, name TEXT)") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE contig(id INTEGER PRIMARY KEY not NULL, name TEXT, strain_id INTEGER, FOREIGN KEY(strain_id) REFERENCES strain(sid))") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE results(id INTEGER PRIMARY KEY not NULL, type TEXT, value TEXT, number INTEGER, start_bp TEXT, end_bp TEXT, locus_id INTEGER, contig_id INTEGER, FOREIGN KEY(locus_id) REFERENCES locus(id),FOREIGN KEY(contig_id) REFERENCES contig(id))") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE locus(id INTEGER PRIMARY KEY not NULL, name TEXT, sequence TEXT, pan TEXT)") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE allele(id INTEGER PRIMARY KEY not NULL, sequence TEXT, copy INTEGER, locus_id INTEGER, contig_id INTEGER, FOREIGN KEY(locus_id) REFERENCES locus(id), FOREIGN KEY(contig_id) REFERENCES contig(id))") or $self->logger->logdie($dbh->errstr);
	
	$dbh->disconnect();
	$self->_sqlString->{'results'}=[];
	$self->_sqlString->{'locus'}=[];
	$self->_sqlString->{'allele'}=[];
	$self->_sqlString->{'strain'}=[];
	$self->_sqlString->{'contig'}=[];
	$self->logger->info("SQL DB initialized");
}

=head2 _sqliteDb

Used for temp file store / retrieval of the core / accessory creation.

=cut
sub _sqliteDb{
	my $self = shift;
	$self->{'__sqliteDb'} = shift // return $self->{'__sqliteDb'};
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

sub _sqlSelectNumber{
	my $self=shift;
	$self->{'__sqlSelectNumber'}=shift // return $self->{'__sqlSelectNumber'};
}

sub _sqlString{
	my $self=shift;
	$self->{'__sqlString'}=shift // return $self->{'__sqlString'};
}

sub _contigIds{
	my $self=shift;
	$self->{'_contigIds'}=shift // return $self->{'_contigIds'};
}

sub _locusId{
	my $self=shift;
	$self->{'_locusId'}=shift // return $self->{'_locusId'};
}

sub _populateStrainTable{
	my $self=shift;
	my $mfsn=shift;

	my %contigIds;
	my $counter=1;

	#my $dbh = (DBI->connect("dbi:SQLite:dbname=" . $self->settings->baseDirectory . "temp_sql.db","","")) or $self->logdie("Could not connect to SQLite DB");
	$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->settings->baseDirectory . "temp_sql.db","","")) or $self->logdie("Could not create SQLite DB");
	
	my $strainId=0;
	my @names = (sort keys %{$mfsn->sequenceNameHash});

	foreach my $name(@names){
		$strainId++;		

		$self->_insertIntoDb(
			table=>'strain',
			id=>$strainId,
			name=>$name
		);

		#add the no value case
		my $missingContig = 'NA_' . $name;	

		$self->_insertIntoDb(
			table=>'contig',
			strain_id=>$strainId,
			name=>$missingContig
		);
		
		$contigIds{$missingContig}=$counter;
		$counter++;
		
		foreach my $contig(@{$mfsn->sequenceNameHash->{$name}->arrayOfHeaders}){
			$self->_insertIntoDb(
				table=>'contig',
				strain_id=>$strainId,
				name=>$contig
			);
			
			$contigIds{$contig}=$counter;
			$counter++;
		}
	}

	$self->_emptySqlBuffers();
	
	$self->_sqliteDb->disconnect();
	$self->logger->info("Strain table populated");
	return \%contigIds;
}

sub run{
	my $self=shift;
	
	$self->logger->info("Analyzing the pan-genome");
	$self->logger->info("Gathering ordered genome names");
	my ($mfsn,$orderedNames)=$self->_generateOrderedNamesArray();
	$self->_mfsn($mfsn);
	$self->_orderedNames($orderedNames);

	$self->logger->info("Populating the strain table");
	$self->_contigIds($self->_populateStrainTable($self->_mfsn));

	my $forker = Parallel::ForkManager->new($self->settings->numberOfCores);
	my $counter=0;
	#process all XML files
	$self->logger->info("Processing Blast output files");
	foreach my $xml(sort @{$self->xmlFiles}){
		$counter++;
		$forker->start and next;
			$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->settings->baseDirectory . "temp_sql.db","","")) or $self->logdie("Could not create SQLite DB");
			$self->_processBlastXML($xml,$counter);
			#unlink $xml;
			$self->_sqliteDb->disconnect();
		$forker->finish;
	}
	$forker->wait_all_children;
	
	#add entries for query segments that have no Blast hits
	if($self->settings->addMissingQuery){
		$self->logger->debug("queryFile specified as " . $self->settings->queryFile);
		$self->_addQueryWithNoBlastHits();
	}
	
	$self->logger->info("Processing blast output files complete");

	#reopen database handle that has been closed from the forking above
	$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->settings->baseDirectory . "temp_sql.db","","")) or $self->logdie("Could not create SQLite DB");

	$self->_createOutputFile('snp',$self->settings->baseDirectory . 'core_snps.txt');
	$self->_createOutputFile('binary',$self->settings->baseDirectory . 'pan_genome.txt');
	if($self->settings->storeAlleles == 1){
		$self->_createAlleleFiles();
	}
	
	#output the pan-genome with the correct locus IDs
	$self->_createCoreAccessoryGenomes();
	$self->_sqliteDb->disconnect();

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



=head2 _addQueryWithNoBlastHits

If a user-supplied queryFile is used in place of pan-genome generation, the possibility
exists for a "pan" fragment to have a hit in 0 genomes.
This function allows for the outputting of such data.

=cut

sub _addQueryWithNoBlastHits{
	my $self=shift;
	
	$self->logger->info("Adding query fragments with no blast hits");
	#get list of all input query names
	my $queryHashRef = $self->_getPanGenomeHashRef();
	
	$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->settings->baseDirectory . "temp_sql.db","","")) or $self->logdie("Could not create SQLite DB");
	my $sql = "SELECT locus.name FROM locus ORDER BY locus.name ASC";
	my $sth = $self->_sqliteDb->prepare($sql) or $self->logger->logdie($self->_sqliteDb->errstr . "\n$sql");
	$sth->execute() or $self->logger->logdie($self->_sqliteDb->errstr);
	
	my %panHits;
	while(my $row = $sth->fetchrow_arrayref){
		$panHits{$row->[0]}=1;
	}
	
	my $missingCounter=0;
	foreach my $panQuery(sort keys %{$queryHashRef}){	
		unless(defined $panHits{$panQuery}){
			#add result if "pan" fragment is completely absent
			my $id = $self->_getUniqueResultId(777) + $missingCounter;
			
			$self->_insertIntoDb(
				table=>'locus',
				id=>$id,
				name=>$panQuery,
				sequence=>$queryHashRef->{$panQuery},
				pan=>'accessory'
			);
			
			foreach my $name(@{$self->_orderedNames}){
				$self->_insertIntoDb(
					table=>'results',
					type=>'binary',
					contig_id=>$self->_contigIds->{'NA_' . $name},
					locus_id=>$id,
					number=>$id,
					start_bp=>0,
					end_bp=>0,
					value=>'0'
				);
			}
			$missingCounter++;
		}
	}	
	$self->_emptySqlBuffers();
	$self->_sqliteDb->disconnect();
	$self->logger->info("$missingCounter query sequences added with no Blast hits");
}

=head2 _createCoreAccessoryGenomes

Based on the user settings, output a fasta file for each
of the core and accessory genomes.
This is in complement to the panGenome that is output via the
_outputPangenomeLocusIds() sub.

=cut

sub _createCoreAccessoryGenomes{
	my $self=shift;
	
	$self->logger->info("Creating the core and accessory fasta files");
	my $sql=qq{
		SELECT locus.id, locus.name, locus.sequence, locus.pan
		FROM locus
		ORDER BY locus.id ASC
	};
	my $sth = $self->_sqliteDb->prepare($sql) or $self->logger->logdie($self->_sqliteDb->errstr . "\n$sql");
	$sth->execute();
	
	my $coreFH = IO::File->new('>' . $self->settings->baseDirectory . 'coreGenomeFragments.fasta') or die "Could not open file coreGenomeFragments.fasta";
	my $accessoryFH = IO::File->new('>' . $self->settings->baseDirectory . 'accessoryGenomeFragments.fasta') or die "Could not open file accessoryGenomeFragments.fasta";
	my $panFH = IO::File->new('>' . $self->settings->baseDirectory . 'panGenomeFragments.fasta') or die "Could not open file panGenomeFragments.fasta";

	while(my $row = $sth->fetchrow_arrayref){
		my $output = '>lcl|' . $row->[0] . '|' . $row->[1] . "\n" . $row->[2] . "\n";
		if($row->[3] eq 'core'){
			$coreFH->print($output);
		}
		elsif($row->[3] eq 'accessory'){
			$accessoryFH->print($output);
		}
		else{
			$self->logger->logdie("Unknown type $row->[3]");
		}
		$panFH->print($output);
	}
	
	$coreFH->close();
	$accessoryFH->close();
	$panFH->close();
}


=head2 _createAlleleFiles

Given each input locus, create a file with the allele for each genome if present.
Note that we use fetchrow_array rather than fethcrow_arrayref, as the same arrayref is returned each time.
This results in row and nextRow being identical, and things fall apart.

=cut

sub _createAlleleFiles{
	my $self=shift;
	
	$self->logger->info("Creating allele files");
	my $sql=qq{
		SELECT strain.name,locus.name,allele.sequence
		FROM allele
		JOIN contig ON allele.contig_id = contig.id
		JOIN strain ON contig.strain_id = strain.id
		JOIN locus ON allele.locus_id = locus.id
		ORDER BY locus.name,strain.name ASC
	};
	
	my $sth = $self->_sqliteDb->prepare($sql);
	$sth->execute();
	
	my $outFH = IO::File->new('>' . $self->settings->baseDirectory . 'locus_alleles.fasta') or die "Could not open file locus_alleles.fasta";
	my @nextRow;
	my @outputBuffer;
	
	my @row = $sth->fetchrow_array;	
	while($row[0]){		
		@nextRow = $sth->fetchrow_array;
		push @outputBuffer,('>' . $row[0] . "\n" . $row[2] . "\n");
		
	    if((!defined $nextRow[0]) || ($row[1] ne $nextRow[1])){   
	    	$outFH->print("Locus ". $row[1] ."\n"); 	
	    	$outFH->print(@outputBuffer);    	
	    	@outputBuffer=();
	    }  
	   
	    @row=@nextRow;
	}
	$outFH->close();
}


=head3 _createOutputFile

Takes output filename and database table name from function call.

=cut

sub _createOutputFile{
	my $self = shift;
	my $type = shift;
	my $outputFile = shift;

	$self->logger->info("Creating output file $outputFile");

	#print all rows from table
	#SELECT * FROM TableA
	#INNER JOIN TableB
	#ON TableA.name = TableB.name
	my $sql;
	if($self->settings->storeAlleles){
		$sql= qq{
			SELECT results.locus_id,locus.name,strain.name,results.value,results.start_bp,results.end_bp,contig.name
		}
	}
	else{
		$sql = qq{
			SELECT results.locus_id,strain.name,results.value,results.start_bp,results.end_bp,contig.name
		};
	}		
	$sql .=qq{
		FROM results
		JOIN locus ON results.locus_id = locus.id
		JOIN contig ON results.contig_id = contig.id
		JOIN strain ON contig.strain_id = strain.id
		WHERE results.type = '$type'
		ORDER BY locus.name,strain.name,results.start_bp ASC
	};

	my $sth = $self->_sqliteDb->prepare($sql) or $self->logger->logdie($self->_sqliteDb->errstr . "\n$sql");
	$sth->execute() or $self->logger->logdie($self->_sqliteDb->errstr);

	my $outFH = IO::File->new('>' . $outputFile) or $self->logger->logdie("Could not create $outputFile");
	#print header for output file
	if($self->settings->storeAlleles){
		$outFH->print("Locus Id\tLocus Name\tGenome\tAllele\tStart bp\tEnd bp\tContig\n");
	}
	else{
		$outFH->print("Locus Id\tGenome\tAllele\tStart bp\tEnd bp\tContig\n");
	}	
	
	while(my $row = $sth->fetchrow_arrayref){
	    $outFH->print(join("\t",@{$row}) . "\n");
	}
	$outFH->close();	
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


sub _getNamesOfGenomesForThisResult{
	my $self = shift;
	my $result = shift;
	
	my %genomeNames;
	foreach my $res(@{$result}){
		my $sn = Modules::Fasta::SequenceName->new($res->[0]);
		$genomeNames{$sn->name}=1;
	}
	return \%genomeNames;
}

sub _processBlastXML {
	my $self = shift;
	my $blastFile = shift;
	my $counter=shift;

	$self->logger->info("Processing Blast output file $blastFile, counter: $counter");
	#this should guarantee a unique number for each result of every Panseq run on the same machine
	#allows up to 1000 SNPs per result
	$counter = $self->_getUniqueResultId($counter);
	#$self->logger->debug("Gather unique id $counter");
	my $blastResult = Modules::Alignment::BlastResults->new($blastFile,$self->settings);
	$self->logger->debug("About to get the first result");
	my $totalResults=0;
	
	while(my $result = $blastResult->getNextResult){
		$totalResults++;
		$self->logger->debug("Result $totalResults processed from BlastResult.pm");
		$self->logger->debug("Result array contains " . scalar(@{$result}) . " elements");
		my $allNames = $self->_getNamesOfGenomesForThisResult($result);
		
		my $numberOfResults = scalar keys %{$allNames};
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
	
		$self->_insertIntoDb(
			table=>'locus',
			id=>$counter,
			name=>$result->[0]->[1],
			sequence=>$result->[0]->[11],
			pan=>$coreOrAccessory
		);
		
		#generate a MSA of all strains that contain sequence
		#if the locus is core, we will send this MSA to the SNPFinder
		my $msaHash = $self->_getHashOfFastaAlignment(
			$self->_getMsa($result,$counter)
		);
		
		foreach my $name(@{$self->_orderedNames}){			
			if(defined $allNames->{$name}){
				next;
			}			
			my $contigId = $self->_contigIds->{'NA_' . $name};
			$self->_insertIntoDb(
						table=>'results',
						type=>'binary',
						contig_id=>$contigId,
						locus_id=>$counter,
						number=>$counter,
						start_bp=>0,
						end_bp=>0,
						value=>'0'
					);
		}
		
		foreach my $res(@{$result}){
			my $contigId = $self->_contigIds->{$res->[0]};
			#if we generate a pan-genome, alleles need to be stored by setting storeAlleles 1 in the config file
			if($self->settings->storeAlleles){
				$self->_insertIntoDb(
					table=>'allele',
					locus_id=>$counter,
					contig_id=>$contigId,
					sequence=>$msaHash->{$res->[0]}
				);
			}					
										
			$self->_insertIntoDb(
				table=>'results',
				type=>'binary',
				contig_id=>$contigId,
				locus_id=>$counter,
				number=>$counter,
				start_bp=>$res->[2],
				end_bp=>$res->[3],
				value=>1
			);				
		}	
	
		if($coreOrAccessory eq 'core'){
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
				$self->logger->debug("contig_id: " . $self->_contigIds->{$cResult->{'contig'}} . "\n"
				. "number: " . $cResult->{'locusId'} . "\n"
				. "locus_id: $counter\n"
				. "start_bp: " .  $cResult->{'startBp'} . "\n"
				. "end_bp: " . $cResult->{'startBp'} . "\n"
				. "value: " . $cResult->{'value'} . "\n");
							
				$self->_insertIntoDb(
					table=>'results',
					type=>'snp',
					contig_id=>$self->_contigIds->{$cResult->{'contig'}},
					number=>$cResult->{'locusId'},
					locus_id=>$counter,
					start_bp=>$cResult->{'startBp'},
					end_bp=>$cResult->{'startBp'},
					value=>$cResult->{'value'}
				);
			}
		}
	}
	$self->logger->info("Total results: $totalResults");
	$self->_emptySqlBuffers();
}


=head2 _getMsa

Take in a Modules::Alignment::BlastResult and generate a MSA of all hit sequences.

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
	foreach my $res(@{$result}){
		$tempInFH->print('>' . $res->[0] . "\n" . $res->[10] . "\n");
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
	
	my %results;
	my $header;
	my $sequence;
	my $alignmentLength;
	
	my %copyNumber=();
	
	foreach my $line(@{$alignedFastaSequences}){
		$line =~ s/\R//g;
		
		if($line =~ /^>(.+)/){		
			$line =~ s/>//;			
			$header= $line;
			
			if(defined $copyNumber{$header}){
				$copyNumber{$header}++;
			}
			else{
				$copyNumber{$header}=1;
			}
			
			$results{$header}->{'copy'}=$copyNumber{$header};
		}
		else{
			if(defined $results{$header}){
				$results{$header}->{'sequence'} .= $line;	
			}
			else{
				$results{$header}->{'sequence'}=$line;
			}					
		}
	}
	return (\%results);
}


=head 2 _emptySqlBuffers

We store the SQL statements to do 500 INSERT at a time.
If any are left in the buffer, add them to the DB.

=cut

sub _emptySqlBuffers{
	my $self=shift;
	
	#the reverse sort is because we need strain to come before contig,
	#as contig uses the strain ID as a foreign key
	foreach my $table(reverse sort keys %{$self->_sqlString}){
		if(defined $self->_sqlString->{$table}->[0]){
			#$self->logger->debug("$table defined, emptying buffers");
			my $sqlString = join('',@{$self->_sqlString->{$table}});
			$self->_sqliteDb->do($sqlString) or $self->logger->logdie("Could not perform SQL DO " . $self->_sqliteDb->errstr . "\n$sqlString");
			#actually empty the buffer, rahter than just print it out
			$self->_sqlString->{$table}=[];
		}
		else{
			#$self->logger->debug("$table not defined, not emptying");
		}
	}	
}

=head2 _insertIntoDb

Take in a table name and any key/value pairs for DB entry.
This is a generalized function that will not check anything beyond the fact
that a parameter named 'table' exists.
Stores all the entries in @keys, without 'table'.

=cut


sub _insertIntoDb{
	my $self = shift;
	my %params = @_;
	
	my @keys;
	foreach my $key(sort keys %params){
		if($key eq 'table'){
			next;
		}
		push @keys, $key;
	}
	my $table = $params{'table'} // $self->logger->logdie("table required in _insertIntoDb");
	#$self->logger->debug("\nTable: $table");
	#taken from http://stackoverflow.com/questions/1609637/is-it-possible-to-insert-multiple-rows-at-a-time-in-an-sqlite-database
	# INSERT INTO 'tablename'
 	# SELECT 'data1' AS 'column1', 'data2' AS 'column2'
	# UNION SELECT 'data3', 'data4'
	# UNION SELECT 'data5', 'data6'
	# UNION SELECT 'data7', 'data8'
	# used UNION ALL due to performance increase (see comments of above linked thread)
	# note that SQLite has default of SQLITE_MAX_COMPOUND_SELECT=500, so we need to account for this

	# Customer
	# ==================
	# Customer_ID | Name

	# Order
	# ==============================
	# Order_ID | Customer_ID | Price
	# insert into "order" (customer_id, price) values \
	# ((select customer_id from customer where name = 'John'), 12.34);

	#http://www.daniweb.com/software-development/perl/threads/339979/assign-query-result-to-variable
	#The following query should return only one row containing one integer.
	# $sth=$dbh->prepare("SELECT 5 as lucky_number;") ||
	# die "Prepare failed: $DBI::errstr\n";
	# $sth->execute() ||
	# die "Couldn't execute query: $DBI::errstr\n";
	# my @record = $sth->fetchrow_array; #Assign result to array variable

	my $sql=[];
	if(defined $self->_sqlString->{$table}->[0]){
		$sql = $self->_sqlString->{$table};
		my $currentSql = "UNION ALL SELECT ";
		
		my $counter=0;
		foreach my $key(@keys){
			#$self->logger->debug("In SQL, key: $key, value: $params{$key}");
			if($counter > 0){
				$currentSql .= ", ";
			}
			
			$currentSql .= "'$params{$key}'";
			$counter++;
		}
		$currentSql .="\n";
		
		#$self->logger->debug("Adding to $table currentSql\n$currentSql");
		push @{$sql},$currentSql;
	}
	else{		
		my $currentSql = "INSERT INTO '$table' (";
		my $counter=0;
		foreach my $key(@keys){
			if($counter > 0){
				$currentSql .= ", ";
			}
			$currentSql .= "$key";
			$counter++;
		}
		$currentSql .= ")\n";
		
		$currentSql .= "SELECT ";
		
		$counter=0;
		foreach my $key(@keys){
			if($counter > 0){
				$currentSql .= ", ";
			}
			
			$currentSql .= "'$params{$key}' AS '$key'";
			$counter++;
		}
		$currentSql .="\n";
		
		#$self->logger->debug("Adding to $table initialSql\n$currentSql");
		push @{$sql}, $currentSql;
		#push @{$sql}, qq{INSERT INTO '$table' (value,start_bp,contig_id,locus_id) SELECT '$value' AS 'value', '$startBp' AS 'start_bp', '$contigId' AS 'contig_id', '$locusId' AS 'locus_id'};
	}
	
	if(scalar(@{$sql})==500){
		my $sqlString = join('',@{$sql});
		#$self->logger->debug("performing DO with $sqlString");
		$self->_sqliteDb->do($sqlString) or $self->logger->logdie("Could not perform SQL DO " . $self->_sqliteDb->errstr . "\n$sqlString");
		$self->_sqlString->{$table}=[];
	}
	else{
		#$self->logger->debug("Adding SQL string to array for table $table");
		$self->_sqlString->{$table}=$sql;
	}
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
	foreach my $res(@{$result}){
#		my $sn = Modules::Fasta::SequenceName->new($res->[0]);
		$startBpHash{$res->[0]}=$res->[4];
	}
	
	#add SNP information to the return
	my $snpDetective = Modules::Alignment::SNPFinder->new(
		 'orderedNames'=>$self->_orderedNames,
		 'alignedFastaHash'=>$msaHash,
		 'startBpHashRef'=>\%startBpHash,
		 'resultNumber'=>$resultNumber,
	 );	
	 my $snpDataArrayRef = $snpDetective->findSNPs();
	 #this returns undef if there are no SNPs
	 return $snpDataArrayRef;	
}

1;





