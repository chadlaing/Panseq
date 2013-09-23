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
    foreach my $key(keys %params){
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
		$self->panGenomeOutputFile($self->outputDirectory . 'pan_genome.txt');
	}

	$self->_currentResult(0);
}

sub _initDb{
	my $self = shift;
	
	$self->logger->info("Initializing SQL DB");
	#define SQLite db
	my $dbh = (DBI->connect("dbi:SQLite:dbname=" . $self->outputDirectory . "temp_sql.db","","")) or $self->logdie("Could not connect to SQLite DB");
	
	$dbh->do("DROP TABLE IF EXISTS results");
	$dbh->do("DROP TABLE IF EXISTS strain");
	$dbh->do("DROP TABLE IF EXISTS contig");
	$dbh->do("DROP TABLE IF EXISTS locus");
	$dbh->do("DROP TABLE IF EXISTS allele");

	$dbh->do("CREATE TABLE strain(id INTEGER PRIMARY KEY not NULL, name TEXT)") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE contig(id INTEGER PRIMARY KEY not NULL, name TEXT, strain_id INTEGER, FOREIGN KEY(strain_id) REFERENCES strain(id))") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE results(id INTEGER PRIMARY KEY not NULL, type TEXT, value TEXT, start_bp TEXT, end_bp TEXT, locus_id INTEGER, contig_id INTEGER, FOREIGN KEY(locus_id) REFERENCES locus(id),FOREIGN KEY(contig_id) REFERENCES contig(id))") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE locus(id INTEGER PRIMARY KEY not NULL, name TEXT)") or $self->logger->logdie($dbh->errstr);
	$dbh->do("CREATE TABLE allele(id INTEGER PRIMARY KEY not NULL, sequence TEXT, locus_id INTEGER, contig_id INTEGER, FOREIGN KEY(locus_id) REFERENCES locus(id), FOREIGN KEY(contig_id) REFERENCES contig(id))") or $self->logger->logdie($dbh->errstr);
	
	$dbh->disconnect();
	$self->_sqlString->{'results'}=[];
	$self->_sqlString->{'locus'}=[];
	$self->_sqlString->{'allele'}=[];
}

=head2 _sqliteDb

Used for temp file store / retrieval of the core / accessory creation.

=cut
sub _sqliteDb{
	my $self = shift;
	$self->{'__sqliteDb'} = shift // return $self->{'__sqliteDb'};
}


=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

=head2 useSuppliedLabels

If the user supplied their own input file, instead of using a generated pan-genome,
we want to keep the supplied loci names.
If this is set to 1, it uses supplied labels. Default is 0.

=cut


sub useSuppliedLabels{
	my $self=shift;
	$self->{'_useSuppliedLabels'}=shift // return $self->{'_useSuppliedLabels'};
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

	my $dbh = (DBI->connect("dbi:SQLite:dbname=" . $self->outputDirectory . "temp_sql.db","","")) or $self->logdie("Could not connect to SQLite DB");
	
	#add case for missing values
	my $sql = qq{INSERT INTO strain(name) VALUES('')};
	$dbh->do($sql);

	foreach my $name(sort keys %{$mfsn->sequenceNameHash}){
		my $sql = qq{INSERT INTO strain(name) VALUES("$name")};
		$dbh->do($sql);

		#add the no value case
		my $missingContig = 'NA_' . $name;
		my $missingSql = qq{INSERT INTO contig(strain_id,name) VALUES((SELECT MAX(id) FROM strain),"$missingContig")};
		$dbh->do($missingSql);
		$contigIds{$missingContig}=$counter;
		$counter++;

		foreach my $contig(@{$mfsn->sequenceNameHash->{$name}->arrayOfHeaders}){
			my $contigSql = qq{INSERT INTO contig(strain_id,name) VALUES((SELECT MAX(id) FROM strain),"$contig")};
			$dbh->do($contigSql);
			$contigIds{$contig}=$counter;
			$counter++;
		}		
	}
	$dbh->disconnect();
	return \%contigIds;
}

sub run{
	my $self=shift;
	
	my ($mfsn,$orderedNames)=$self->_generateOrderedNamesArray();
	$self->_mfsn($mfsn);
	$self->_orderedNames($orderedNames);

	$self->_contigIds($self->_populateStrainTable($self->_mfsn));

	my $forker = Parallel::ForkManager->new($self->numberOfCores);
	my $counter=0;
	#process all XML files
	foreach my $xml(@{$self->xmlFiles}){
		$counter++;
		$forker->start and next;
			$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->outputDirectory . "temp_sql.db","","")) or $self->logdie("Could not create SQLite DB");
			$self->_processBlastXML($xml,$counter);
			unlink $xml;
			$self->_sqliteDb->disconnect();
		$forker->finish;
	}
	$forker->wait_all_children;
	$self->logger->info("Processing blast output files complete");

	for my $count(1..2){
		$forker->start and next;
		#reopen database handle that has been closed from the forking above
		$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->outputDirectory . "temp_sql.db","","")) or $self->logdie("Could not create SQLite DB");
		if($count==1){
			$self->_createOutputFile('snp',$self->outputDirectory . 'core_snps.txt');
		}
		elsif($count==2){
			$self->_createOutputFile('binary',$self->outputDirectory . 'pan_genome.txt');
			if($self->useSuppliedLabels){
				$self->_createAlleleFiles();
			}
		}
		else{
			$self->logger->logconfess("Count is $count, but should not be greater than 2");
		}
		$self->_sqliteDb->disconnect();
		$forker->finish;
	}
	$forker->wait_all_children;
	$self->logger->info("Pan-genome generation complete");
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
	
	my $outFH = IO::File->new('>' . $self->outputDirectory . 'locus_alleles.fasta') or die "Could not open file locus_alleles.fasta";
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
	my $sql = qq{
		SELECT results.locus_id,locus.name,strain.name,results.value,results.start_bp,results.end_bp,contig.name
		FROM results
		JOIN locus ON results.locus_id = locus.id
		JOIN contig ON results.contig_id = contig.id
		JOIN strain ON contig.strain_id = strain.id
		WHERE results.type = '$type'
		ORDER BY locus.name,results.locus_id,strain.name,results.start_bp,results.end_bp,contig.name ASC
	};

	my $sth = $self->_sqliteDb->prepare($sql) or $self->logger->logdie($self->_sqliteDb->errstr . "\n$sql");
	$sth->execute() or $self->logger->logdie($self->_sqliteDb->errstr);

	my $outFH = IO::File->new('>' . $outputFile) or $self->logger->logdie("Could not create $outputFile");
	#print header for output file
	$outFH->print("Locus Id\tLocus Name\tGenome\tAllele\tStart bp\tEnd bp\tContig\n");
	
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


sub _processBlastXML {
	my $self = shift;
	my $blastFile = shift;
	my $counter=shift;

	$self->logger->info("Processing Blast output file $blastFile, counter: $counter");
	$counter *=1000000000;

	my $blastResult = Modules::Alignment::BlastResults->new($blastFile,$self->percentIdentityCutoff);
	
	while(my $result = $blastResult->getNextResult){
		my @names = keys %{$result};
		$counter++;		
		
		$self->_insertIntoDb(
			table=>'locus',
			id=>$counter,
			name=>$result->{$names[0]}->[1]
		);
		
		#generate a MSA of all strains that contain sequence
		#if the locus is core, we will send this MSA to the SNPFinder
		my $msaHash = $self->_getHashOfFastaAlignment(
			$self->_getMsa($result,$counter)
		);
		
		foreach my $name(@{$self->_orderedNames}){							
			if(defined $result->{$name}){
					#currently, we only store the alleles if using an input file
					#if we generate a pan-genome, no alleles are stored
					if($self->useSuppliedLabels){
						$self->_insertIntoDb(
							table=>'allele',
							locus_id=>$counter,
							contig_id=>$self->_contigIds->{$result->{$name}->[0]},
							sequence=>$msaHash->{$name}->{'sequence'}
						);
					}					
										
					$self->_insertIntoDb(
						table=>'results',
						type=>'binary',
						contig_id=>$self->_contigIds->{$result->{$name}->[0]},
						locus_id=>$counter,
						start_bp=>$result->{$name}->[2],
						end_bp=>$result->{$name}->[3],
						value=>1
					);
				}
				else{
					$self->_insertIntoDb(
						table=>'results',
						type=>'binary',
						contig_id=>$self->_contigIds->{'NA_' . $name},
						locus_id=>$counter,
						start_bp=>0,
						end_bp=>0,
						value=>'0'
					);
				}		
			}	
	
		if(scalar(keys %{$result}) >= $self->coreGenomeThreshold){
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
				$self->_insertIntoDb(
					table=>'results',
					type=>'snp',
					contig_id=>$self->_contigIds->{$cResult->{'contig'}},
					locus_id=>$cResult->{'locusId'},
					start_bp=>$cResult->{'startBp'},
					end_bp=>$cResult->{'startBp'},
					value=>$cResult->{'value'}
				);
			}
		}
	}
	$self->logger->info("Total results: $counter");
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
	my $tempInFile = $self->outputDirectory . 'muscleTemp_in' . $resultNumber;
	my $tempInFH = IO::File->new('>'. $tempInFile) or die "$!";
	my $tempOutFile = $self->outputDirectory . 'muscleTemp_out' . $resultNumber;
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
	foreach my $hit(sort keys %{$result}){
		$tempInFH->print('>' . $result->{$hit}->[0] . "\n" . $result->{$hit}->[10] . "\n");
	}
	
	my $systemLine = $self->muscleExecutable . ' -in ' . $tempInFile . ' -out ' . $tempOutFile . ' -maxiters 3 -quiet';
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
$hashRef->{'name'}->{'fasta'}="fasta header"
$hashRef->{'name'}->{'sequence'}="DNA sequence"
the {'fasta'} key contains the fasta header,
and the {'sequence'} key contains the DNA sequence.

=cut

sub _getHashOfFastaAlignment{
	my $self = shift;
	my $alignedFastaSequences=shift;
	
	my $results={};
	my $name;
	my $alignmentLength;

	foreach my $line(@{$alignedFastaSequences}){
		$line =~ s/\R//g;
		
		if($line =~ /^>(.+)/){		
			$line =~ s/>//;
			my $sn = Modules::Fasta::SequenceName->new($line);
			$name = $sn->name;
			$self->logger->debug("New name: $name");
			#we don't need or want the '>'; all names are stored without the fasta header signifier
			$results->{$name}->{'fasta'}=$line;
		}
		else{
			if(defined $results->{$name}->{'sequence'}){
			#	$self->logger->debug("Sequence already present, adding to it:\n$line");
				$results->{$name}->{'sequence'}=$results->{$name}->{'sequence'} . $line;	
			}
			else{
				#$self->logger->debug("New sequence, adding:\n$line");
				$results->{$name}->{'sequence'}=$line;
			}					
		}
	}
	return ($results);
}


=head 2 _emptySqlBuffers

We store the SQL statements to do 500 INSERT at a time.
If any are left in the buffer, add them to the DB.

=cut

sub _emptySqlBuffers{
	my $self=shift;
	
	foreach my $table(keys %{$self->_sqlString}){
		if(defined $self->_sqlString->{$table}->[0]){
			my $sqlString = join('',@{$self->_sqlString->{$table}});
			$self->_sqliteDb->do($sqlString) or $self->logger->logdie("Could not perform SQL DO " . $self->_sqliteDb->errstr);
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
			if($counter > 0){
				$currentSql .= ", ";
			}
			
			$currentSql .= "'$params{$key}'";
			$counter++;
		}
		$currentSql .="\n";
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
		push @{$sql}, $currentSql;
		#push @{$sql}, qq{INSERT INTO '$table' (value,start_bp,contig_id,locus_id) SELECT '$value' AS 'value', '$startBp' AS 'start_bp', '$contigId' AS 'contig_id', '$locusId' AS 'locus_id'};
	}
	
	if(scalar(@{$sql})==500){
		my $sqlString = join('',@{$sql});
		#$self->logger->info("$sqlString");
		$self->_sqliteDb->do($sqlString) or $self->logger->logdie("Could not perform SQL DO " . $self->_sqliteDb->errstr . "\n$sqlString");
		$self->_sqlString->{$table}=[];
	}
	else{
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
	foreach my $hit(sort keys %{$result}){
		$startBpHash{$result->{$hit}->[0]}=$result->{$hit}->[4];
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





