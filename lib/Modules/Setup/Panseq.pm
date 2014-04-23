#!/usr/bin/env perl

=pod

=head1 NAME

Modules::Setup::Panseq 

=head1 SYNOPSIS

=head1 DESCRIPTION

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

package Modules::Setup::Panseq;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use File::Basename;
use Modules::Setup::PanseqFiles;
use Modules::NovelRegion::NovelIterator;
use Modules::Alignment::MakeBlastDB;
use Modules::Fasta::SegmentMaker;
use Modules::Fasta::FastaFileSplitter;
use Modules::PanGenome::PanGenome;
use Modules::LociSelector::LociSelector;
use Parallel::ForkManager;
use Tie::Log4perl;
use Log::Log4perl;
use File::Path qw/make_path remove_tree/;
use Carp;
use File::Copy;
use Archive::Zip;
use Role::Tiny::With;
use DBI;

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

    my $settingsObj = shift;
   
    #get and parse the configuration file for Panseq
	$self->settings($settingsObj);
	$self->_createDirectories();

	#logging
    $self->logger(Log::Log4perl->get_logger()); 
	$self->logger->info("Logger initialized in Modules::Setup::Panseq");  
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

=head3 settings

All Panseq settings from config file.

=cut

sub settings{
	my $self=shift;
	$self->{'_settings'}=shift // return $self->{'_settings'};
}

sub run{
	my $self=shift;

	#if mode is set to loci, launch loci_finder, rather than Panseq
	if($self->settings->runMode eq 'loci'){
		$self->_launchLociFinder();
	}
	else{
		$self->_launchPanseq();
	}
	#$self->_cleanUp();
	$self->_createZipFile();
}

=head2 _launchLociFinder

Runs the loci_finder.pl script.

=cut

sub _launchLociFinder{
	my $self = shift;
	
	my $files = Modules::Setup::PanseqFiles->new(
		'queryDirectory'=>$self->settings->queryDirectory,
		'referenceDirectory'=>$self->settings->referenceDirectory // undef
	);
	
	my $numberOfLoci;
	if($self->settings->numberOfLoci ==0){
		$numberOfLoci='best';
	}
	else{
		$numberOfLoci = $self->settings->numberOfLoci;
	}

	
	if(defined $files->queryFileNames->[0]){
		my $lociLine = "perl $FindBin::Bin/loci_selector.pl " 
		. $files->queryFileNames->[0]
		. ' '
		. $numberOfLoci . ' > ' . $self->settings->baseDirectory . "selectedLoci.txt";
		system($lociLine);
		
	}
	else{
		$self->logger->logdie("No input file for _launchLociFinder");
	}
}


=head2 _launchPanseq

Runs either pan-genome or novel region finding, based on the $self->settings->runMode option.
If a queryFile is specified, it will be used as the 'panGenomeFile' instead of the $files->singleQueryFile,
and it will also skip the novel region gather. However, $files->singleQueryFile will be used as the
database file that the 'panGenomeFile' is compared to.
If $self->settings->queryFile is defined, it is used as the pan-genome file regardless of other settings.

=cut

sub _launchPanseq{
	my $self = shift;

	#get the query/reference files to use in the comparison
	my $novelIterator;
	my $files = Modules::Setup::PanseqFiles->new(
		'queryDirectory'=>$self->settings->queryDirectory,
		'referenceDirectory'=>$self->settings->referenceDirectory // undef
	);

	if(defined $self->settings->queryFile){
		#sanitize the input file for proper fasta format
		#use the sanitized file for future work
		#with Roles::CombineFilesIntoSingleFile->_combineAndSanitizeFastaFiles

		#requires ($arrayRef, outputFileName)
		my $sanitizedFileName = $self->settings->baseDirectory . 'sanitized_queryFile.fasta';
		$self->_combineAndSanitizeFastaFiles(
			[$self->settings->queryFile],
			$sanitizedFileName
		);
		
		$novelIterator=Modules::NovelRegion::NovelIterator->new(
			'panGenomeFile'=>$sanitizedFileName,
			'queryFile'=>$files->singleQueryFile($self->settings->baseDirectory . 'singleQueryFile.fasta'),
			'settings'=>$self->settings
		); 
	}
	else{
		#get novel regions
		$novelIterator = Modules::NovelRegion::NovelIterator->new(
			'queryFile'=>$files->singleQueryFile($self->settings->baseDirectory . 'singleQueryFile.fasta'),
			'referenceFile'=>$files->singleReferenceFile($self->settings->baseDirectory . 'singleReferenceFile.fasta'),
			'novelRegionsFile'=>$self->settings->baseDirectory . 'novelRegions.fasta',
			'settings'=>$self->settings
		);
		$novelIterator->run();
	}

	$self->logger->info("Panseq mode set as " . $self->settings->runMode);
	#perform pan-genomic analyses
	if(defined $self->settings->runMode && $self->settings->runMode eq 'pan'){
		$self->_performPanGenomeAnalyses($files,$novelIterator);
		#with File::Copy
		move($novelIterator->panGenomeFile,$self->settings->baseDirectory . 'panGenome.fasta');
		$self->_createTreeFiles();
	}else{
		#with File::Copy
		move($novelIterator->panGenomeFile,$self->settings->baseDirectory . 'novelRegions.fasta');
	}
}


=head3

Given the $self->settings->baseDirectory,
test to see if it exists.
If it does, stop program with a message.
If it does not, create directory, as well as directory/logs for logging output.

=cut

sub _createDirectories{
	my $self=shift;

	if(-d $self->settings->baseDirectory){
		print "\nThe directory you have specified for program output:\n\t" . $self->settings->baseDirectory . "\n"
			. "already exists.\n";
		
		print "overwrite: " . $self->settings->overwrite . "\n";
		
		if($self->settings->overwrite){
			print "Overwrite set to true. Deleting directory " . $self->settings->baseDirectory . "\n";
			remove_tree($self->settings->baseDirectory);
		}
		else{		
			print "Please specify a new baseDirectory\n";
			exit(1);
		}
	}
	make_path($self->settings->baseDirectory);
	
	unless(-d $self->settings->baseDirectory){
		print STDERR "Unable to create directory " . $self->settings->baseDirectory . "\n";
		exit(1);
	}
}


=head3 _cleanUp

Removes any files left-over from the run.
Uses Roles::CombineFilesIntoSingleFile for _getFileNamesFromDirectory

=cut

sub _cleanUp{
	my $self=shift;
	
	my $fileNamesRef = $self->_getFileNamesFromDirectory($self->settings->baseDirectory);

	foreach my $file(@{$fileNamesRef}){
		if(
			($file =~ m/_db(\.n|temp)/) || 
			($file =~ m/_sequenceSplitterDB/) || 
			($file =~ m/singleQueryFile\.fasta/) ||
			($file =~ m/\.delta$/) ||
			($file =~ m/(accessory|core|muscle|nucmer)Temp/) ||
			($file =~ m/\.xml/) ||
			($file =~ m/ReferenceFile/) ||
			($file =~ m/_withRefDirectory_temp/) ||
			($file =~ m/lastNovelRegionsFile/) ||
			($file =~ m/uniqueNovelRegions/) ||
			($file =~ m/temp_sql\.db/) ||
			($file =~ m/\.temp/)
		){
			unlink $file;
		}
	}
}

=head3 _createTreeFiles

Create both a core-snp and pan-genome +/- file for use in phylogenetic analyses.
If more than one processor available, make both at the same time.

=cut

sub _createTreeFiles{
	my $self = shift;

	my $forker= Parallel::ForkManager->new($self->settings->numberOfCores);

	for my $num(1..2){
		$forker->start and next;
			if($num==1){
				$self->_createTree('snp');
			}
			elsif($num==2){
				$self->_createTree('binary');
			}
			else{
				$self->logger->logconfess("num value is $num, should not exceed 2");
			}
		$forker->finish();
	}
	$forker->wait_all_children();
}

sub _createTree{
	my $self=shift;
	my $type = shift;
	
	$self->logger->info("Creating tree $type");
	#define SQLite db
	my $dbh = (DBI->connect("dbi:SQLite:dbname=" . $self->settings->baseDirectory . "temp_sql.db","","")) or $self->logger->logdie("Could not connect to SQLite DB");
	
	my $sql;
	$self->logger->info("Name or id: " . $self->settings->nameOrId);

	if($self->settings->nameOrId eq 'name'){
		$sql = "SELECT strain.name,results.value, locus.name, results.number";
	}
	elsif($self->settings->nameOrId eq 'id'){
		$sql = "SELECT strain.name,results.value, locus.id, results.number";
	}
	else{
		$self->logger->logdie("Incorrect setting " . $self->settings->nameOrId . " in nameOrId. Expected name or id");
	}
	
	$sql .= qq{
		FROM results
		JOIN contig ON results.contig_id = contig.id
		JOIN strain ON contig.strain_id = strain.id
		JOIN locus ON results.locus_id = locus.id
		WHERE results.type = '$type'
		ORDER BY locus.name, results.number ASC
	};

	my $sth = $dbh->prepare($sql);
	my $didSucceed = $sth->execute();

	unless($didSucceed){
		$self->logger->fatal($self->logger->logdie("sth->execute fail" . $dbh->errstr . "\n$sql"));
	}
	
	my $tableFH = IO::File->new('>' . $self->settings->baseDirectory . $type . '_table.txt') or $self->logger->logdie("$!");
	my %results;
	my %loci;
	my $locus;
	my @genomeOrder;
	
	my @row = $sth->fetchrow_array;
	my $counter=0;
	while(defined $row[0]){
		$counter++;
	  	my @nextRow = $sth->fetchrow_array;
	  	$loci{$row[0]}=$row[1];
	  	
	    push @{$results{$row[0]}},$row[1];

	    if((!defined $nextRow[0]) || ($row[3] ne $nextRow[3])){
	    	unless(defined $genomeOrder[0]){
	    		@genomeOrder = sort keys %loci;
	    		$tableFH->print("\t" . join("\t",@genomeOrder) . "\n");
	    	}
	    	$tableFH->print($row[2]);
	    	foreach my $genome(@genomeOrder){
	    		$tableFH->print("\t" . $loci{$genome});
	    	}
	    	$tableFH->print("\n");
	    }
	    @row=@nextRow;  
	}
	
	
	$dbh->disconnect();
	$tableFH->close();	

	my $nameConversion = $self->_printPhylipFile($type,\%results);
	
	if($type eq 'binary'){
		$self->_printConversionInformation($nameConversion);
	}	
}

sub _printPhylipFile{
	my $self=shift;
	my $table = shift;
	my $results = shift;

	my $outFH = IO::File->new('>' . $self->settings->baseDirectory . $table . '.phylip') or die "$!";

	my $counter=1;
	my %nameConversion;
	foreach my $genome(sort keys %{$results}){
		$nameConversion{$counter}=$genome;

		if($counter==1){
			$outFH->print(scalar(keys %{$results}) . "\t" . scalar(@{$results->{$genome}}) . "\n");
		}

		$outFH->print($counter . $self->_numberOfSpacesToAdd($counter) . join('',@{$results->{$genome}}) . "\n");
		$counter++;
	}
	$outFH->close();
	return \%nameConversion;
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
	my $hashRef =shift;

	my $conversionFH = IO::File->new('>' . $self->settings->baseDirectory . 'phylip_name_conversion.txt') or die "$!";
	$conversionFH->print(
		'Number' . "\t" . 'Name' . "\n"
	);

	foreach my $number(sort keys %{$hashRef}){
		$conversionFH->print($number . "\t" . $hashRef->{$number} . "\n");
	}

	$conversionFH->close();	
}


=head2 _performPanGenomeAnalyses

If a fragmentation size is set, fragment the pan-genome.
Assign either the original or fragmented pan-genome to $panGenomeFile.
The query directory has previously been combined into $files->singleQueryFile,
which is used as the database to screen the pan-genome against. If providing a single
queryFile to represent the panGenome, it should be outside of the queryDirectory, so that it
is not included in the database.
The list of XML files is then processed and the panAnalyzer (Modules::PanGenome::PanGenome)
object is returned.

=cut

sub _performPanGenomeAnalyses{
	my $self=shift;
	my $files=shift;
	my $novelIterator=shift;

	#fragmentationSize defaults to 0 if no fragmentation is to be done
	my $panGenomeFile = $novelIterator->panGenomeFile;
	my $segmenter;

	if($self->settings->fragmentationSize > 0){
		$segmenter = Modules::Fasta::SegmentMaker->new(
			'inputFile'=>$panGenomeFile,
			'outputFile'=>$self->settings->baseDirectory . 'pangenome_fragments.fasta',
			'segmentSize'=>$self->settings->fragmentationSize
		);
		$segmenter->segmentTheSequence;
		$panGenomeFile = $segmenter->outputFile;
	}	

	my $dbCreator = Modules::Alignment::MakeBlastDB->new(
		'dbtype'=>'nucl',
		'blastDirectory'=>$self->settings->blastDirectory,
		'in'=>$files->singleQueryFile,		
		'out'=>$self->settings->baseDirectory . 'panseq_db',
		'title'=>'panseq_db',
		'logfile'=>'>>'.$self->settings->baseDirectory . 'Master.log'	
	);
	$dbCreator->run();
	
	my $splitter= Modules::Fasta::FastaFileSplitter->new(
		inputFile=>$panGenomeFile,
		databaseFile=>$self->settings->baseDirectory . 'fasta_splitter.temp',
		numberOfSplits=>$self->settings->numberOfCores
	);
	$splitter->splitFastaFile();

	my @blastFiles;
	my $forker = Parallel::ForkManager->new($self->settings->numberOfCores);

	foreach my $splitFile(@{$splitter->arrayOfSplitFiles}){
		my $blastOutFile = $splitFile . '_blast.out';
		push @blastFiles, $blastOutFile;
		$forker->start and next;
			my $blastLine = $self->settings->blastDirectory . 'blastn -query ' . $splitFile . ' -out ' . $blastOutFile
				. ' -db ' . $self->settings->baseDirectory . $dbCreator->title 
				. ' -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen pident length sseq qseq"' 
				. ' -evalue 0.001 -word_size 20 -num_threads 1'
				. ' -max_target_seqs 100000';
			$self->logger->info("Running blast with the following: $blastLine");
			system($blastLine);
			#unlink $splitFile;
		$forker->finish;
	}
	$forker->wait_all_children();
	
	#do the pan-genome analysis	
	my $panAnalyzer = Modules::PanGenome::PanGenome->new(
		xmlFiles=>\@blastFiles,
		queryFile=>$files->singleQueryFile,
		panGenome=>$panGenomeFile,
		settings=>$self->settings
	);
	$panAnalyzer->run();
	
	if(defined $segmenter && -e $segmenter->outputFile){
		unlink $segmenter->outputFile;
	}
	
	#add the functional assignment of pan-genome loci based of blastx of the genbank NR database
	#only if nrDatabase is defined
	if(defined $self->settings->nrDatabase){
		$self->_assignFunction($panGenomeFile);
	}
	
	return 1;
}


=head2 _assignFunction

Using blastx and the NR database, assign a function to each query locus.

=cut

sub _assignFunction{
	my $self=shift;
	my $queryFile=shift;	
	
	my $blastLine = $self->settings->blastDirectory . 'blastx -query ' . $queryFile . ' -out ' . $self->settings->baseDirectory . 'blastx_nr.out'
		. ' -db ' . $self->settings->nrDatabase
		. ' -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen pident length sseq qseq stitle"' 
		. ' -evalue 0.001 -word_size 3 -num_threads ' . $self->settings->numberOfCores
		. ' -max_target_seqs 1 -gilist ' . "$FindBin::Bin/bacteria_gi_list";
	$self->logger->info("Running blastx with the following: $blastLine");
	system($blastLine);
}


sub _createZipFile{
    my $self=shift;
    
    $self->logger->info("Creating zip file");
    #object to add all files the user receives
    my $zipper = Archive::Zip->new();
    
    #this is used in creating the batch file for the program run. called 'baseDirectory' within the Panseq scripts.
    my $outputDir =  $self->settings->baseDirectory;

    #get all file names from directory (to avoid dealing with logging, not using FileManipulation.pm)
    opendir( DIRECTORY, $outputDir ) or die "cannot open directory $outputDir $!\n";
    my @dir = readdir DIRECTORY;
    closedir DIRECTORY;
    foreach my $fileName (@dir) {
            #dont add directories or . and .. files
            next if substr( $fileName, 0, 1 ) eq '.';                

            #usage is: addFile( $fileName [, $newName ] ) from Archive::Zip manual
            my $fullName = $outputDir . $fileName;
            $zipper->addFile($fullName, $fileName);
    }
    my $status = $zipper->writeToFileNamed($outputDir . 'panseq_results.zip');
    
    if($status == 0){
            #everything ok
            return 1;
    }
    else{
            #something went wrong, return error
            return 0;
    }
}

1;
