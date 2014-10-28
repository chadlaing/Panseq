#!/usr/bin/env perl

=pod

=head1 NAME

Modules::NovelRegion::NovelRegionFinder - Finds novel regions for a genome or group of genomes when compared to a genome
or group of genomes.

=head1 SYNOPSIS


	use Modules::NovelRegion::NovelRegionFinder;
	
	my $nrf = Modules::NovelRegion::NovelRegionFinder->new(
		'mode'=>$self->settings->novelRegionFinderMode,
		'coordsFile'=>$nucmer->coordsFile,
		'queryFile'=>$nucmer->queryFile,
		'minimumNovelRegionSize'=>$self->settings->minimumNovelRegionSize
	);
	$nrf->findNovelRegions();
	$nrf->printNovelRegions($retriever, $novelOutputFile);

=head1 DESCRIPTION

Finds and stores the novel regions for a genome or group of genomes in the following format: 
$hashRef->id->,1..5,8..134,678..45999 etc.
There are three modes of operation:
'no_duplicates' - regions present in one or more query sequences but no reference sequences.
'unique' - regions present only in one of the query strains and no others.
'common_to_all' - regions present in all of the query strains and no reference strains.

=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: 

=head1 AUTHOR

Chad Laing (chadlaing@gmail.com)

=head2 Methods

=cut

package Modules::NovelRegion::NovelRegionFinder;

#includes
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Bio::SeqIO;
use Carp;
use Log::Log4perl;
use Role::Tiny::With;

with 'Roles::GetNewIdStartEnd';

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

    $self->logger->debug("Logger initialized in Modules::NovelRegion::NovelRegionFinder");  

    my %params = @_;

    #on object construction set all parameters
    foreach my $key(sort keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::NovelRegion::NovelRegionFinder");
		}
	}	


	#initialize anonymous hashes and arrays
	$self->_queryFastaHeadersHash({});
	$self->_comparisonHash({});

	#set default parameters
	unless(defined $self->adjacentJoiningSize){
		$self->adjacentJoiningSize(10);
	}
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

=head2 settings

An object containing all of the settings from the user.

=cut

sub settings{
	my $self=shift;
	$self->{'_settings'}=shift // return $self->{'_settings'};
}


=head3 coordsFile

The .coords file produced from the nucmer run.
It is used to get coordinates to determine novel regions.

=cut

sub coordsFile{
	my $self=shift;
	$self->{'_coordsFile'}=shift // return $self->{'_coordsFile'};
}



=head3 queryFile

The file containing all combined query sequences.
File is in mult-fasta foramt.

=cut

sub queryFile{
	my $self=shift;
	$self->{'_queryFile'}=shift // return $self->{'_queryFile'};
}

=head3 _queryFastaHeadersHash

A hash in the form of {headerName}=lengthOfSequence, for every fasta sequence in the multi-fasta file.

=cut

sub _queryFastaHeadersHash{
	my $self=shift;
	$self->{'__queryFastaHeadersHash'}=shift // return $self->{'__queryFastaHeadersHash'};
}


=head3 _comparisonHash

Stores the alignment between fasta sequences in the form {query}->{reference}=,start..end,start..end

=cut

sub _comparisonHash{
	my $self=shift;
	$self->{'__comparisonHash'}=shift // return $self->{'__comparisonHash'};
}


=head3 _novelRegionsHashRef

Stores the coordinates of all the novel regions given the query file and the mummer coords file.

=cut

sub novelRegionsHashRef{
	my $self=shift;
	$self->{'__novelRegionsHashRef'}=shift // return $self->{'__novelRegionsHashRef'};
}


=head3 mode

Sets whether a no_duplicates, unique, or common_to_all run is run.

=cut

sub mode{
	my $self=shift;
	$self->{'_mode'} = shift // return $self->{'_mode'};
}

=head3 adjacentJoiningSize

On a single contig, if two regions are novel and separated by a space, they will be joined provided the space
is less than or equal to this parameter. If not set, it defaults to 10 bp.

=cut

sub adjacentJoiningSize{
	my $self=shift;
	$self->{'_adjacentJoiningSize'} = shift // return $self->{'_adjacentJoiningSize'};
}


sub findNovelRegions {
	my ($self) = shift;


	#construct hash of all fasta headers from all query files in _queryFastaHeadersHash
	#populate the _quryNameObjectHash
	$self->_processQueryFiles();

	$self->_buildHashOfAllComparisons();

	#get the novel regions and set the class variable 
	$self->novelRegionsHashRef($self->_dispatcher);
}

sub _processQueryFiles{
	my $self=shift;

	$self->logger->debug("Gathering fasta sequence information from:" . $self->queryFile);

	my $inFH = Bio::SeqIO->new(-file=>'<'.$self->queryFile, -format=>'fasta') or $self->logger->logdie("Cannot open " . $self->queryFile . "$!");

	while(my $seq = $inFH->next_seq){
		my $header = $seq->id . $seq->desc;
		$self->_queryFastaHeadersHash->{$header}=$seq->length();
	}
	$inFH->close();
}

sub _buildHashOfAllComparisons {
	my ($self) = shift;

	$self->logger->debug("Building hash of all comparisons");

	my $inputFH = IO::File->new( '<' . $self->coordsFile) or $self->logger->logdie("Cannot open " .$self->coordsFile . " $!");

	while (my $line = $inputFH->getline) {

		$self->_updateComparisonHash($line);
	}
	$inputFH->close();
	$self->_printComparisonHash();
}


=head3 printNovelRegions

	->printNovelRegions(
		$retrieverObject,$outputFile
	);
	
novelRegionsHashRef is in the form:
$hashRef->id->,1..5,8..134,678..45999 etc.
The id is the fasta id indexed in the database and the string of locations ,X..Y represent a single subsequence to return.
This method calls the extractRegion method and prints the results greater than or equal to 'cutoffSize' to outputFile.

=cut

sub printNovelRegions{
	my $self=shift;
	my $retriever = shift // $self->logger->logconfess("A retriever object is required in printNovelRegions");
	my $outputFile = shift // $self->logger->logconfess("An output filename is required in printNovelRegions");

	my $outFH=IO::File->new('>' . $outputFile) or $self->logger->logdie("Could not open $outputFile");
		
	foreach my $id(sort keys %{$self->novelRegionsHashRef}){
		my $coordString = $self->novelRegionsHashRef->{$id};
			
		while($coordString =~ /\,(\d+)\.\.(\d+)/gc){				
			my $start =$1;
			my $end =$2;
			my $length=$end-$start+1;
			
			$self->logger->debug("Length: $length cutoff: " . $self->settings->minimumNovelRegionSize);
			next unless $length >= $self->settings->minimumNovelRegionSize;	

			#uses Roles::GetNewIdStartEnd to implement _getNewIdStartEnd
			my($relId, $relStart, $relEnd) = $self->_getNewIdStartEnd($id,$start,$end);
			my $relLength = ($relEnd - $relStart +1);

			$self->logger->debug("original: start:$start, end:$end, length: $length\nnew    : start:$relStart, end:$relEnd, length:$relLength" );
			$outFH->print('>' . $relId . '_(' . $relStart . '..' . $relEnd . ')'. "\n" . $retriever->extractRegion($id,$start,$end) . "\n");	
		}# end while
	}# end foreach
	$outFH->close();	
}



=head3 _updateComparisonHash

Parses the current line of the nucmer coords file.
Adds to any current hits, creates new ones if two fasta headers have never been matched.
All the alignments are based on the fasta header name, not the SequenceName

=cut

sub _updateComparisonHash {
	my ($self) = shift;
	my $line = shift;

	#remove end of line characters
	$line =~ s/\R//g;

	my @la = split(/\t/,$line);

	#we need to only deal with actual results, not headers
	#return and do nothing if not a result
	unless(scalar(@la) == 11){
		return;
	}

	# my %hit=(
	# 	'rStart'=>$la[0],
	# 	'rEnd'=>$la[1],
	# 	'qStart'=>$la[2],
	# 	'qEnd'=>$la[3],
	# 	'rHitLength'=>$la[4],
	# 	'qHitLength'=>$la[5],
	# 	'percentId'=>$la[6],
	# 	'rSeqLength'=>$la[7],
	# 	'qSeqLength'=>$la[8],
	# 	'rName'=>$la[9],
	# 	'qName'=>$la[10]
	# );

	#start,end,key,value
	$self->_addToComparisonHash($la[2],$la[3],$la[10],$la[9]);

	# if($self->settings->novelRegionFinderMode eq 'unique'){
	# 	$self->logger->debug("Adding REF:SEQ for unique");
	# 	$self->_addToComparisonHash($la[0],$la[1],$la[9],$la[10]);
	# }
}


=head2 _addToComparisonHash

We want to be able to use query or reference as the keys to the hash.
This allows the construction of unique regions much easier and efficient (avoiding N^2).

=cut

sub _addToComparisonHash{
	my $self = shift;
	my $start =shift;
	my $end = shift;
	my $keyName = shift;
	my $valueName = shift;

	if($start > $end){
		$self->logger->debug("Swapping $start with $end");

		#perl's handy swap number function
		($start,$end) =($end,$start);
	}

	my $alignmentCoords;
	if(defined $self->_comparisonHash->{$keyName}->{$valueName}){
		$self->logger->debug("Alignment coords previously exist for $keyName:$valueName");
		$alignmentCoords = $self->_comparisonHash->{$keyName}->{$valueName};
	}

	if(defined $alignmentCoords){
		$self->logger->debug("Before update: $start\.\.$end -- $alignmentCoords");
		$self->_comparisonHash->{$keyName}->{$valueName}=$self->_updateAlignmentCoords($alignmentCoords, $start, $end);
	}
	else{
		$alignmentCoords = ',' . $start . '..' . $end;

		$self->logger->debug("New alignment: $keyName:$valueName $alignmentCoords");

		$self->_comparisonHash->{$keyName}->{$valueName}=$alignmentCoords;
	}
}



sub _updateAlignmentCoords {
	my $self = shift;
	my $alignmentCoords = shift;
	my $qStart          = shift;    #current values 'q', checking whether they can modify the start/end of any coords
	my $qEnd            = shift;

	$self->logger->debug("qStart: $qStart qEnd: $qEnd" );

	my $newStart;
	my $newEnd;
	my $prevStart;
	my $prevEnd;
	my $lastTime  = 0;
	my $duplicate = 0;
	my $oldString;
	my $newString;

	while ( $alignmentCoords =~ /\,(\d+)\.\.(\d+)/gc ) {
		$prevStart = $1;
		$prevEnd   = $2;

		$self->logger->debug("prevStart: $prevStart prevEnd: $prevEnd" );

		#check for a q start that is less than or between the current values
		if ( $qStart > $prevEnd ) {
			next;
		}

		#if both qStart and qEnd are within the current set, no need to add the duplicate
		if ( ( $qStart >= $prevStart ) && ( $qEnd <= $prevEnd ) ) {
			$duplicate = 1;
			last;
		}

		#initial values
		$newStart = $prevStart;
		$newEnd   = $prevEnd;

		#extend the beginning range if needed
		if ( ( $qStart < $prevStart ) && ( $qEnd >= $prevStart ) ) {
			$newStart = $qStart;
			$lastTime = 1;
		}

		#extend the end range if needed
		if ( ( $qEnd > $prevEnd ) && ( $qStart <= $prevEnd ) ) {
			$newEnd   = $qEnd;
			$lastTime = 1;
		}

		if ($lastTime) {
			$oldString = ',' . $prevStart . '..' . $prevEnd;
			last;
		}
	}    #end while

	#double check something has changed, then update
	if ( defined $oldString ) {
		$newString = ',' . $newStart . '..' . $newEnd;

		$self->logger->debug("Replacing $oldString with ,$newStart\.\.$newEnd");

		$alignmentCoords =~ s/\Q$oldString\E/$newString/;
	}
	elsif ( !$duplicate ) {

		#this means nothing was modified, but the new values need to be added
		$newString = ',' . $qStart . '..' . $qEnd;

		$self->logger->debug("Adding new match: $newString");
		$alignmentCoords .= $newString;
	}

	#otherwise, duplicate value, nothing added
	$self->logger->debug("Updated alignmentCoords: $alignmentCoords" );
	return $alignmentCoords;
}

sub _getNoDuplicates {
	my ($self) = shift;

	my %nonDuplicatedCoords;	

	#algorithm
	#A vs refs
	#B vs A + refs
	#C vs. A + B + refs
	#...
	
	foreach my $query ( sort keys %{$self->_queryFastaHeadersHash}) {
		$self->logger->debug("query: $query");
		my $qName = $self->settings->getGenomeNameFromContig($query);
		
		if($qName eq '_ignore'){
			next;
		}
		
		#in case there are no hits, we want the query contig to be present in the output
		unless(defined $self->_comparisonHash->{$query}){
			$nonDuplicatedCoords{$query} = ',0..0';
		}
		
		foreach my $ref (sort keys %{ $self->_comparisonHash->{$query} } ) {
			
			my $rName = $self->settings->getGenomeNameFromContig($ref);
			#the only reason the same sequence should be in the query and ref file
			if($qName eq $rName){
				$self->logger->debug("query $query is the same");
				unless(defined $nonDuplicatedCoords{$query}){
					$nonDuplicatedCoords{$query} = ',0..0';
				}
				next;
			}
			
			$nonDuplicatedCoords{$query} .= $self->_comparisonHash->{$query}->{$ref};
			$self->logger->debug("query: $query ref: $ref non_dup_coords: $nonDuplicatedCoords{$query}");
		}
	}
	$self->logger->debug("Sorting and compiling novel regions");
	my $sortedHashRef = $self->_sortAndCompileCoordsHash( \%nonDuplicatedCoords );

	$self->logger->debug("Getting negative image of coords");
	return $self->_getNegativeImageCoords($sortedHashRef);
}

sub _sortAndCompileCoordsHash {
	my ($self) = shift;
	my $hashRef = shift;

	foreach my $query (sort keys %{$hashRef} ) {
		$self->logger->debug("Sorting and compiling $query" );
		$hashRef->{$query} = $self->_sortAndCompileCoordString($hashRef->{$query} );
	}
	return $hashRef;
}


sub _sortAndCompileCoordString {
	my ($self) = shift;

	if (@_) {
		my $coords = shift;
		my %coordsHash;
		my $newString;

		$self->logger->debug("Coords: $coords");
		#create hash with first coord as key
		while ( $coords =~ /(\,(\d+)\.\.\d+)/gc ) {

			#do the sorting and compiling
			my $coord    = $1;
			my $startKey = $2;

			#account for duplicate start coords --> otherwise problems
			while ( defined $coordsHash{$startKey} ) {
				$startKey++;
			}
			$coordsHash{$startKey} = $coord;
		}

		#sort hash, reconstruct string
		foreach my $key ( sort { $a <=> $b } keys %coordsHash ) {
			if ( defined $newString ) {
				my $start;
				my $end;
				my $sortCoord;
				if ( $coordsHash{$key} =~ /(\,(\d+)\.\.(\d+))/ ) {
					$sortCoord = $1;
					$start     = $2;
					$end       = $3;

					$self->logger->debug("\tSort and compile: Before update: Key: $key start: $start end: $end" );

				}
				else {
					$self->logger->logconfess($coordsHash{$key} . ' did not match!');
				}

				$newString = $self->_updateAlignmentCoords( $newString, $start, $end );
				$self->logger->debug("Sort and compile: After update: Key: $key start: $start end: $end" );
			}
			else {
				$newString = $coordsHash{$key};
			}
		}
		$self->logger->debug("Sort and compile: New string: $newString" );
		return $newString;
	}
	else {
		$self->logconfess("nothing sent to sortAndCompileCoordString");
	}
}

sub _getNegativeImageCoords {
	my ($self) = shift;

	my $hashRef = shift;
	my %negativeHash;

	#hash of $hash{$query}=<coords in order>
	foreach my $key (sort keys %{$hashRef} ) {
		my $missingCoords;
		my $prevEnd = 0;
		
		#get the negative image coords
		my $nextKey = $hashRef->{$key};
		
		while ( $nextKey =~ /\,(\d+)\.\.(\d+)/gc ) {
			my $start = $1;
			my $end   = $2;

			$self->logger->debug("start: $start, end: $end");

			if ( ( $start > ($prevEnd + 1) ) && ( ( $start + 1 ) < $end )) {

				$missingCoords .= ',' . ( $prevEnd + 1 ) . '..' . ( $start - 1 );

				$self->logger->debug("NEG:\tpassed checks: missingCoords: $missingCoords prevEnd:$prevEnd start:$start");
			}
			$prevEnd = $end;
		}

		#check the last end
		if ( ( $self->_queryFastaHeadersHash->{$key} > $prevEnd )) {

			#if the whole thing is missing, set start to 1, rather than 0
			# if($prevEnd ==0){
			# 	$prevEnd=1;
			# }
			$missingCoords .= ',' . ($prevEnd+1) . '..' . ( $self->_queryFastaHeadersHash->{$key} );
		}

		$self->logger->debug("NEG: original of $key: $nextKey");
		if ( defined $missingCoords ) {
			$self->logger->debug("NEG: negative: $missingCoords");
		}
		else {
			$self->logger->debug("NEG: no missing coords");
		}

		$negativeHash{$key} = $missingCoords if defined $missingCoords;
	}
	return \%negativeHash;

}


sub _dispatcher {
	my ($self) = shift;

	if ( defined $self->settings->novelRegionFinderMode ) {

		my $novelCoordsHashRef;
		if ( $self->settings->novelRegionFinderMode eq 'no_duplicates' ) {
			$self->logger->debug("Getting non-redundant query regions");
			$novelCoordsHashRef = $self->_getNoDuplicates();
		}
		#elsif ( $self->novelRegionFinderMode eq 'common_to_all' ) {
		# 	$novelCoordsHashRef = $self->getCommonToAll();
		# }
		elsif ( $self->settings->novelRegionFinderMode eq 'unique' ) {
			$self->logger->info("Getting unique regions of query strains");
		 	$novelCoordsHashRef = $self->_getNoDuplicates();
		}
		else {
		 	$self->logger->logconfess("incorrect type sent to getNovelRegions");
		}		
		return $self->_joinAdjacentNovelRegionsHash($novelCoordsHashRef);
	}
	else {
		$self->logger->logconfess("comparisonType not defined in getNovelRegions");
	}
}

sub _joinAdjacentNovelRegionsHash {
	my $self = shift;

	if ( scalar(@_) == 1 ) {
		my $hashRef = shift;

		my %joinedHash;

		foreach my $seq (sort keys %{$hashRef} ) {
			unless ( defined $hashRef->{$seq} ) {
				$self->logger->debug("$seq has no novel coords");
				next;
			}
			$joinedHash{$seq} = $self->_joinAdjacentNovelRegionsString( $hashRef->{$seq} );
		}
		return \%joinedHash;
	}
	else {
		$self->logger->logconfess("incorrect number of arguments sent to joinAdjacentNovelRegions\n");
	}
}

sub _joinAdjacentNovelRegionsString {
	my ($self) = shift;

	my $coordsString = shift // $self->logger->logconfess("UNDEF coordsString sent to _joinAdjacentNovelRegionsString");

	$self->logger->debug("coordsString in joinAdjacentNovelRegionsString is: $coordsString");
	return $coordsString if ( $self->adjacentJoiningSize == '0' );

	my $prevStart;
	my $prevEnd;
	my $start;
	my $end;
	my $prevCoords;
	my $currCoords;
	my $newCoordsString = $coordsString;

	my $count = 0;
	while ( $coordsString =~ /(\,(\d+)\.\.(\d+))/gc ) {

		$currCoords = $1;
		$start      = $2;
		$end        = $3;

		next unless defined $prevStart;

		if ( ( $start - $prevEnd ) < $self->adjacentJoiningSize ) {
			my $oldCoords = $prevCoords . $currCoords;
			my $newCoords = ',' . $prevStart . '..' . $end;
			$newCoordsString =~ s/$oldCoords/$newCoords/;
		}
	}
	continue {
		$prevStart  = $start;
		$prevEnd    = $end;
		$prevCoords = $currCoords;
	}
	return $newCoordsString;
}



sub _printComparisonHash{
	my $self=shift;
	
	foreach my $key(sort keys %{$self->_comparisonHash}){
		foreach my $key2 (sort keys %{$self->_comparisonHash->{$key}}){
			$self->logger->debug($key . ' : ' . $key2 . ' : ' . $self->_comparisonHash->{$key}->{$key2});
		}
	}
}




1;
