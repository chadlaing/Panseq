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
	$self->_locusId({});
	$self->_initDb();

	#default values
	unless(defined $self->panGenomeOutputFile){
		$self->panGenomeOutputFile($self->outputDirectory . 'pan_genome.txt');
	}

	unless(defined $self->coreSnpsOutputFile){
		$self->coreSnpsOutputFile($self->outputDirectory . 'core_snps.txt');
	}

	$self->_currentResult(0);
}

sub _initDb{
	my $self = shift;

	#define SQLite db
	my $dbh = (DBI->connect("dbi:SQLite:dbname=" . $self->outputDirectory . "temp_sql.db","","")) or $self->logdie("Could not connect to SQLite DB");
	
	$dbh->do("DROP TABLE IF EXISTS snp");
	$dbh->do("DROP TABLE IF EXISTS binary");
	$dbh->do("DROP TABLE IF EXISTS strain");
	$dbh->do("DROP TABLE IF EXISTS contig");

	$dbh->do("CREATE TABLE strain(id INTEGER PRIMARY KEY, name TEXT)");
	$dbh->do("CREATE TABLE contig(id INTEGER PRIMARY KEY, name TEXT, strain_id INTEGER, FOREIGN KEY(strain_id) REFERENCES strain(id))");
	$dbh->do("CREATE TABLE snp(id INTEGER PRIMARY KEY, value TEXT, start_bp TEXT, locus_id INTEGER, contig_id INTEGER, FOREIGN KEY(contig_id) REFERENCES contig(id))");
	$dbh->do("CREATE TABLE binary(id INTEGER PRIMARY KEY, value TEXT, start_bp TEXT,locus_id INTEGER, contig_id INTEGER, FOREIGN KEY(contig_id) REFERENCES contig(id))");

	$dbh->disconnect();

	$self->_sqlString->{'binary'}=[];
	$self->_sqlString->{'snp'}=[];
	$self->_locusId->{'snp'}=1;
	$self->_locusId->{'binary'}=1;
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
			$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->outputDirectory . "temp_sql.db","","")) or $self->logdie("Could not vreate SQLite DB");
			$self->_processBlastXML($xml,$counter);
			unlink $xml;
			$self->_sqliteDb->disconnect();
		$forker->finish;
	}
	$forker->wait_all_children;
	$self->logger->info("Procssing XML files complete");

	for my $count(1..2){
		$forker->start and next;
		#reopen database handle that has been closed from the forking above
		$self->_sqliteDb(DBI->connect("dbi:SQLite:dbname=" . $self->outputDirectory . "temp_sql.db","","")) or $self->logdie("Could not vreate SQLite DB");
		if($count==1){
			$self->_createOutputFile('snp',$self->outputDirectory . 'core_snps.txt');
		}
		elsif($count==2){
			$self->_createOutputFile('binary',$self->outputDirectory . 'pan_genome.txt');
		}
		else{
			$self->logger->logconfess("Count is $count, but should not be greater than 2");
		}	
		$forker->finish;
	}
	$forker->wait_all_children;
	$self->logger->info("Pan-genome generation complete");
}

=head3 _createOutputFile

Takes output filename and database table name from function call.

=cut

sub _createOutputFile{
	my $self = shift;
	my $table = shift;
	my $outputFile = shift;

	$self->logger->info("Creating output file $outputFile");

	#print all rows from table
	#SELECT * FROM TableA
	#INNER JOIN TableB
	#ON TableA.name = TableB.name
	my $sql = qq{
		SELECT $table.locus_id,strain.name,$table.value,$table.start_bp,contig.name 
		FROM $table
		JOIN contig ON $table.contig_id = contig.id
		JOIN strain ON contig.strain_id = strain.id 
	};
	#my $sql = qq{SELECT locus_id,value,start_bp,contig_id FROM $table};
	my $sth = $self->_sqliteDb->prepare($sql);
	$sth->execute();

	my $outFH = IO::File->new('>' . $outputFile) or $self->logger->logdie("Could not create $outputFile");
	#print header for output file
	$outFH->print("Locus Id\tGenome\tAllele\tStart bp\tContig\n");
	
	while(my $row = $sth->fetchrow_arrayref){
	    $outFH->print(join("\t",@{$row}) . "\n");
	}
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
	
	return ($multiFastaSN,[sort keys %{$multiFastaSN->sequenceNameHash}])
	#order the hash for a consistent order of names
}


sub _processBlastXML {
	my $self = shift;
	my $blastXMLFile = shift;
	my $counter=shift;

	$self->logger->info("Processing Blast XML file $blastXMLFile, counter: $counter");
	$counter *=1000000000;

	my $xmlFH = IO::File->new( '<' . $blastXMLFile ) or die "$!";
	my $xmler = Modules::Alignment::BlastResultFactory->new($xmlFH);

	while ( my $result = $xmler->nextResult) {
		$counter++;
		my $type;
		if($self->coreGenomeThreshold == 0){
			$type = 'accessory';
		}
		else{
			$type = $self->_getCoreAccessoryType($result);
		}
		$self->_processResult($result, $type,$counter);
		
	}	
	$xmlFH->close();

	#process any remaining sql that didn't fill up the 500 select buffer
	if(defined $self->_sqlString->{'binary'}->[0]){
		my $sqlString = join('',@{$self->_sqlString->{'binary'}});
		$self->_sqliteDb->do($sqlString) or $self->logger->logdie("$!");
	}
	
	if(defined $self->_sqlString->{'snp'}->[0]){
		my $sqlString = join('',@{$self->_sqlString->{'snp'}});
		$self->_sqliteDb->do($sqlString) or $self->logger->logdie("$!");
	}
}

sub _insertIntoDb{
	my $self = shift;
	my $table = shift;
	my $contigId = shift;
	my $locusId = shift;
	my $startBp = shift;
	my $value = shift;

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
		push @{$sql}, qq{ UNION ALL SELECT '$value','$startBp', '$contigId','$locusId'};
	}
	else{
		push @{$sql}, qq{INSERT INTO '$table' (value,start_bp,contig_id,locus_id) SELECT '$value' AS 'value', '$startBp' AS 'start_bp', '$contigId' AS 'contig_id', '$locusId' AS 'locus_id'};
	}
	
	if(scalar(@{$sql})==500){
		my $sqlString = join('',@{$sql});
		$self->_sqliteDb->do($sqlString) or $self->logger->logdie("$!");
		$self->_sqlString->{$table}=[];
	}
	else{
		$self->_sqlString->{$table}=$sql;
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
			$self->logger->debug("coverage: $sequenceCoverage" . " cutoff: " . $self->percentIdentityCutoff);
		}

		if ( $numberOverSequenceCutoff >= $self->coreGenomeThreshold ) {
			$returnType = 'core';
		}
		else {
			$returnType = 'accessory';
		}
	}
	$self->logger->debug("$returnType number over cutoff: $numberOverSequenceCutoff");
	$self->logger->debug("Returning $returnType");
	return $returnType;
}


sub _getSequenceCoverage {
	my $self = shift;

	my $hit         = shift // $self->logger->logconfess('hit is required in Modules::PanGenome::_getSequenceCoverage');
	my $queryLength = shift // $self->logger->logconfess('queryLength is required in Modules::PanGenome::_getSequenceCoverage');

	my $sequenceCoverage = ( $hit->hsp_align_len - ( $hit->hsp_align_len - $hit->hsp_identity ) ) / $queryLength * 100;
	return $sequenceCoverage;
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
	my $counter=shift;

	my @fastaArray;
	my %startBpHash;
	$counter *=1000000;

	if ( defined $result->hitHash ) {
	
		#each $hit is the SequenceName::Name of the hit
		foreach my $name (@{$self->_orderedNames} ) {
			my $hitObj;
			
			if(defined $result->hitHash->{$name}){
				$self->logger->info("$name found in $counter hitHash");
				$hitObj = $result->hitHash->{$name};
				$hitObj->setSequence();    #this ensures start <  end bp

				#table, contigId,locusId,startBp,value
				$self->_insertIntoDb(
					'binary',
					$self->_contigIds->{$hitObj->hit_def},
					$counter,
					$hitObj->hsp_hit_from,
					$self->_getAccessoryValue(
						$self->_getSequenceCoverage( $hitObj, $result->query_len ),
						$hitObj
					)
				);				
			}
			else{
				$self->logger->info("$name not found in $counter hitHash");
				#table, contigId,locusId,startBp,value
				$self->_insertIntoDb(
					'binary',
					$self->_contigIds->{"NA_$name"},
					$counter,
					'0',
					'0'
				);
				next;
			}				
			
			#create a fasta-file format array if core
			if($type eq 'core'){		
 				push @fastaArray, ( '>' . $hitObj->hit_def . "\n" . $hitObj->hsp_hseq . "\n");
 				$startBpHash{ $hitObj->hit_def } = $hitObj->hsp_hit_from;
			}
		}    #end of foreach
		if($type eq 'core'){
			my $coreResultArray = $self->_getCoreResult(\@fastaArray, \%startBpHash,$counter,$self->_locusId->{'snp'});

			foreach my $result(@{$coreResultArray}){
				$self->_insertIntoDb(
					'snp',
					$self->_contigIds->{$result->{'contig'}},
					$result->{'locusId'},
					$result->{'startBp'},
					$result->{'value'}
				);
			}
		}
	}#end of if
	else{
		$self->logger->info("Result :" . $result->name . " has no hits!");
	}
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

sub _getCoreResult {
	my $self = shift;

	my $fastaArrayRef = shift;
	my $startBpHashRef = shift;	
	my $resultNumber=shift;

	#create temp files for muscle
	my $tempInFile = $self->outputDirectory . 'muscleTemp_in' . $resultNumber;
	my $tempInFH = IO::File->new('>'. $tempInFile) or die "$!";
	my $tempOutFile = $self->outputDirectory . 'muscleTemp_out' . $resultNumber;
	my $tempOutFH = IO::File->new('+>' . $tempOutFile) or die "$!";

	$tempInFH->print(@{$fastaArrayRef});	
	my $systemLine = $self->muscleExecutable . ' -in ' . $tempInFile . ' -out ' . $tempOutFile . ' -maxiters 3 -quiet';
	system($systemLine);

	# #close the open FH
	$tempInFH->close();	
	my @alignedFastaSeqs = $tempOutFH->getlines();
	$tempOutFH->close();

	# #delete temp files
	unlink $tempInFile;
	unlink $tempOutFile;

	#add SNP information to the return
	my $snpDetective = Modules::Alignment::SNPFinder->new(
		 'orderedNames'=>$self->_orderedNames,
		 'alignedFastaSequences'=>\@alignedFastaSeqs,
		 'startBpHashRef'=>$startBpHashRef,
		 'resultNumber'=>$resultNumber,
	 );	
	 my $snpDataArrayRef = $snpDetective->findSNPs();
	 #this returns undef if there are no SNPs
	 return $snpDataArrayRef;	
}



1;





