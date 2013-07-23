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
use diagnostics;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Bio::SeqIO;
use Modules::Fasta::SequenceName;
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
    foreach my $key(keys %params){
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
	$self->_queryNameObjectHash({});
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


=head3 coordsFile

The .coords file produced from the nucmer run.
It is used to get coordinates to determine novel regions.

=cut

sub coordsFile{
	my $self=shift;
	$self->{'_coordsFile'}=shift // return $self->{'_coordsFile'};
}


=head3 minimumNovelRegionSize

Specifies the size in bp that a novel region must be greater than or equal to.
Only those novel regions passing this threshold will be out put to the results file.

=cut

sub minimumNovelRegionSize{
	my $self=shift;
	$self->{'_minimumNovelRegionSize'}=shift // return $self->{'_minimumNovelRegionSize'};
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


=head3 _queryNameObjectHash

A hash in the form {NameObject->name}=NameObject, where
NameObject are those from Modules::Fasta::SequenceName.

=cut

sub _queryNameObjectHash{
	my $self=shift;
	$self->{'__queryNameObjectHash'}=shift // return $self->{'__queryNameObjectHash'};
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

	#build hash of all comparisons
	$self->_buildHashOfAllComparisons();

	#get the novel regions and set the class variable 
	$self->novelRegionsHashRef($self->_dispatcher);
}

sub _processQueryFiles{
	my $self=shift;

	$self->logger->info("Gathering fasta sequence information from:" . $self->queryFile);

	my $inFH = Bio::SeqIO->new(-file=>'<'.$self->queryFile, -format=>'fasta') or $self->logger->logdie("Cannot open " . $self->queryFile . "$!");

	while(my $seq = $inFH->next_seq){
		my $header = $seq->id . $seq->desc;
		my $sn = Modules::Fasta::SequenceName->new($header);
		$self->_queryFastaHeadersHash->{$header}=$seq->length();
		$self->_queryNameObjectHash->{$sn->name}=$sn;

		$self->logger->debug("Processing $header into: " . $sn->name . ' length: ' . $seq->length);
	}
	$inFH->close();
}

sub _buildHashOfAllComparisons {
	my ($self) = shift;

	$self->logger->info("Building hash of all comparisons");

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
		
	foreach my $id(keys %{$self->novelRegionsHashRef}){
		my $coordString = $self->novelRegionsHashRef->{$id};
			
		while($coordString =~ /\,(\d+)\.\.(\d+)/gc){				
			my $start =$1;
			my $end =$2;
			my $length=$end-$start+1;
			
			#$self->logger->info("Length: $length cutoff: " . $self->minimumNovelRegionSize);

			next unless $length >= $self->minimumNovelRegionSize;	

			#uses Roles::GetNewIdStartEnd to implement _getNewIdStartEnd
			my($relId, $relStart, $relEnd) = $self->_getNewIdStartEnd($id,$start,$end);
		
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

	my %hit=(
		'rStart'=>$la[0],
		'rEnd'=>$la[1],
		'qStart'=>$la[2],
		'qEnd'=>$la[3],
		'rHitLength'=>$la[4],
		'qHitLength'=>$la[5],
		'percentId'=>$la[6],
		'rSeqLength'=>$la[7],
		'qSeqLength'=>$la[8],
		'rName'=>$la[9],
		'qName'=>$la[10]
	);


	if($hit{'qStart'} > $hit{'qEnd'}){
		$self->logger->debug("Swapping $hit{'qStart'} with $hit{'qEnd'}");

		#perl's handy swap number function
		($hit{'qStart'},$hit{'qEnd'}) =($hit{'qEnd'},$hit{'qStart'});
	}

	my $alignmentCoords;
	if(defined $self->_comparisonHash->{$hit{'qName'}}->{$hit{'rName'}}){
		$self->logger->debug("Alignment coords previously exist for $hit{'qName'}:$hit{'rName'}");
		$alignmentCoords = $self->_comparisonHash->{$hit{'qName'}}->{$hit{'rName'}};
	}

	if(defined $alignmentCoords){
		$self->logger->debug("Before update: $hit{'qStart'}\.\.$hit{'qEnd'} -- $alignmentCoords");
		$self->_comparisonHash->{$hit{'qName'}}->{$hit{'rName'}}=$self->_updateAlignmentCoords($alignmentCoords, $hit{'qStart'}, $hit{'qEnd'});
	}
	else{
		$alignmentCoords = ',' . $hit{'qStart'} . '..' . $hit{'qEnd'};

		$self->logger->debug("New alignment: $hit{'qName'}:$hit{'rName'} $alignmentCoords");

		$self->_comparisonHash->{$hit{'qName'}}->{$hit{'rName'}}=$alignmentCoords;
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
	
	foreach my $query (sort keys %{$self->_queryFastaHeadersHash} ) {
		my $queryName = Modules::Fasta::SequenceName->new($query);
		
		#in case there are no hits, we want the query contig to be present in the output
		unless($self->_comparisonHash->{$query}){
			$nonDuplicatedCoords{$query} = ',0..0';
		}
		
		foreach my $ref ( keys %{ $self->_comparisonHash->{$query} } ) {
			
			my $referenceName = Modules::Fasta::SequenceName->new($ref);
			next if $queryName->name eq $referenceName->name;
			
			$nonDuplicatedCoords{$query} .= $self->_comparisonHash->{$query}->{$ref};
			$self->logger->debug("query: $query ref: $ref non_dup_coords: $nonDuplicatedCoords{$query}\n");
		}
	}
	$self->logger->info("Sorting and compiling novel regions");
	my $sortedHashRef = $self->_sortAndCompileCoordsHash( \%nonDuplicatedCoords );

	$self->logger->info("Getting negative image of coords");
	return $self->_getNegativeImageCoords($sortedHashRef);
}

sub _sortAndCompileCoordsHash {
	my ($self) = shift;
	my $hashRef = shift;

	foreach my $query ( keys %{$hashRef} ) {
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
	foreach my $key ( keys %{$hashRef} ) {
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

	if ( defined $self->mode ) {

		my $novelCoordsHashRef;

		if ( $self->mode eq 'no_duplicates' ) {
			$self->logger->info("Getting non-redundant query regions");
			$novelCoordsHashRef = $self->_getNoDuplicates();
		}
		#elsif ( $self->novelRegionFinderMode eq 'common_to_all' ) {
		# 	$novelCoordsHashRef = $self->getCommonToAll();
		# }
		elsif ( $self->mode eq 'unique' ) {
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


# sub getCommonToAll {
# 	my ($self) = shift;

# 	my %commonCoordsHash;

# 	#use one query sequence as the one from which to extract the sequence data
# 	#common in all allows this
# 	#for each query, check that all other queries have a hit
# 	#using "one" as the w.r.t sequence, get all reference hits for the query
# 	#take negative image

# 	my %oneSeqHash;    #this is a hash so that it can be fed into sortAndCompileCoordsHash

# 	my @queryNames = keys %{ $self->queryNameObjectHash };
# 	my $oneSeq     = FileInteraction::Fasta::SequenceName->new( $queryNames[0] );

# 	$self->logger->debug( "CTA: oneSeq: " . $oneSeq->name );

# 	foreach my $query (sort keys %{ $self->comparisonHash } ) {
# 		$self->logger->debug("CTA: query: $query");

# 		my $querySeqName = FileInteraction::Fasta::SequenceName->new($query);
# 		next unless $querySeqName->name eq $oneSeq->name;

# 		#get the "oneSeq" coords
# 		foreach my $hit ( keys %{ $self->comparisonHash->{$query} } ) {
# 			$self->logger->debug("CTA: hit: $hit");

# 			if ( defined $self->_queryFastaNamesHash->{$hit} ) {
# 				next;
# 				#we only want sequence matches in the reference, for the negative image
# 			}
# 			else {
# 				$self->logger->debug("$hit not defined in queryFastaNamesHash");
# 			}
# 			$oneSeqHash{$query} .= $self->comparisonHash->{$query}->{$hit};

# 			$self->logger->debug( "CTA: oneSeq: " . $oneSeq->name . " oneSeqHash{query}:" . $oneSeqHash{$query} );
# 		}
# 	}

# 	return $self->trimOneSeqToAllQuerySequences( $self->getNegativeImageCoords( $self->sortAndCompileCoordsHash( \%oneSeqHash ) ) );

# }

# sub trimOneSeqToAllQuerySequences {
# 	my ($self) = shift;

# 	if (@_) {
# 		my $hashRef = shift;    #this is the hashRef of all reference regions not found in oneSeq
# 		my %trimmedHash;

# 		$self->logger->debug( "Trim: begin trim, number of hashRef keys is: " . scalar keys %{$hashRef} );

# 		#should only be one sequence name, but possibly multiple keys due to the oneSeq selection in getCommonToAll
# 		foreach my $oneSeqQueryName ( keys %{$hashRef} ) {
# 			my $comparisonHashRef = $self->comparisonHash->{$oneSeqQueryName};

# 			#the following loop stores only query hits in the %coords hash
# 			my %coords;
# 			foreach my $comparisonHit ( keys %{$comparisonHashRef} ) {
# 				$self->logger->debug("Trim: comparisonHit: $comparisonHit");
# 				next unless defined $self->_queryFastaNamesHash->{$comparisonHit};
# 				$self->logger->debug( "Trim: made the cut: $comparisonHit: " . $comparisonHashRef->{$comparisonHit} );
# 				$coords{$comparisonHit} = $comparisonHashRef->{$comparisonHit};
# 			}
# 			$trimmedHash{$oneSeqQueryName} =
# 			  $self->getCommonCoordsWRTinputCoords( $hashRef->{$oneSeqQueryName}, \%coords );    #input: coords as string , hashRefToCheck
# 		}
# 		return \%trimmedHash;
# 	}
# 	else {
# 		confess "nothing sent to trimOneSeqToAllQuerySequences\n";
# 	}
# }

# sub getCommonCoordsWRTinputCoords {
# 	my ($self) = shift;

# 	if ( scalar(@_) == 2 ) {
# 		my $inputCoordString = shift // die "undefined inputCoordString in getCommonCoordsWRTinputCoords\n";
# 		my $hashRef = shift;    #all the query seqs and their hits (no ref hits)
# 		my $newString;

# 		$self->logger->debug("WRT inputCoordString: $inputCoordString");

# 		#algorithm
# 		#get seqName list for each inputCoord (all query sequence matches for the given oneSeq are in the hashRef)
# 		#cycle through seqName list, trimming as it goes

# 		my $seqNameCoordsHashRef = $self->getSeqNameCoordList($hashRef);

# 		my $iteration = 0;

# 		while ( $inputCoordString =~ /(\,\d+\.\.\d+)/gc ) {
# 			my $currentCoords = $1;

# 			$self->logger->debug("WRT currentCoords: $currentCoords");

# 			my $returnedValue = $self->getTrimmedCoords( $currentCoords, $seqNameCoordsHashRef );

# 			if ( defined $returnedValue ) {
# 				$self->logger->debug("WRT afterTrim: $returnedValue");
# 			}
# 			else {
# 				$self->logger->debug("WRT afterTrim: not common to all");
# 				next;
# 			}

# 			$newString .= $returnedValue;
# 			$self->logger->debug("WRT newString: $newString");
# 		}
# 		return $newString;
# 	}
# 	else {
# 		confess "nothing sent to getCommonCoordsWRTinputCoords\n";
# 	}
# }

# sub getTrimmedCoords {
# 	my ($self) = shift;

# 	#This sub creates a string of 0s representing the length of coords sent to the sub 0000000000
# 	#It then increments each position by 1 for every match eg. 0001111100, then 00012222200
# 	#and outputs only the region from the original string that has numbers equal the total number of query sequences

# 	if ( scalar(@_) == 2 ) {
# 		my $coords  = shift;
# 		my $hashRef = shift;    # $hashRef->{seqName}->{seqFastaHeader}=<coords>

# 		#get start/stop coords
# 		my $start;
# 		my $end;

# 		if ( $coords =~ /\,(\d+)\.\.(\d+)/ ) {
# 			$start = $1;
# 			$end   = $2;
# 		}
# 		else {
# 			confess "no coords sent to getTrimmedCoords\n";
# 		}

# 		my $offset = $start - 1;

# 		$self->logger->debug("GTC coords: $coords start: $start end: $end offset:$offset");

# 		my $sequence = '0' . ( '0' x ( $end - $start + 1 ) );    #to account for 1 offset

# 		#check that each seqName has a hit, and where it is
# 		my $seqNameCounter = 0;
# 		my $trigger        = 1;
# 		my $previousName;

# 		foreach my $sName ( keys %{$hashRef} ) {
# 			if ( $trigger == 1 ) {
# 				$seqNameCounter++;
# 				$self->logger->info("Matching $sName");
# 			}
# 			else {
# 				$self->logger->warn("Common sequence not present in getTrimmedCoords for $previousName");
# 				last;
# 			}

# 			$trigger      = 0;
# 			$previousName = $sName;

# 			$self->logger->debug("sName: $sName seqNameCounter: $seqNameCounter");

# 			foreach my $eachMatch ( keys %{ $hashRef->{$sName} } ) {
# 				my $eachCoord = $hashRef->{$sName}->{$eachMatch};

# 				$self->logger->debug("eachMatch $eachMatch");

# 				while ( $eachCoord =~ /\,(\d+)\.\.(\d+)/gc ) {

# 					my $eachStart = $1;
# 					my $eachEnd   = $2;

# 					$self->logger->debug("$eachMatch $eachStart:$eachEnd");

# 					#check to make sure each start/end is within the range
# 					if ( $eachStart > $end ) {
# 						$self->logger->debug("NEXT eachStart greater than end");
# 						next;
# 					}

# 					if ( $eachEnd < $start ) {
# 						$self->logger->debug("NEXT eachEnd less than start");
# 						next;
# 					}
# 					$self->logger->debug("MATCH eachStart: $eachStart eachEnd: $eachEnd");
# 					$eachStart = $start if $eachStart < $start;
# 					$eachEnd   = $end   if $eachEnd > $end;

# 					my $eachLength = $eachEnd - $eachStart + 1;

# 					$self->logger->debug("Before adjustment eachStart: $eachStart start:$start eachEnd: $eachEnd end:$end");

# 					#account for difference between absolute and relative locations
# 					$eachStart = $eachStart - $offset;
# 					$eachEnd   = $eachEnd - $offset;

# 					#get substring, increase count by 1, replace in sequence string
# 					$self->logger->debug("After adjustment eachStart: $eachStart eachEnd: $eachEnd");

# 					my $substring = substr( $sequence, $eachStart, ( $eachEnd - $eachStart + 1 ) );
# 					my $prevNum = $seqNameCounter - 1;

# 					$substring =~ s/$prevNum/$seqNameCounter/g;
# 					substr( $sequence, $eachStart, $eachLength ) = $substring;
# 					$self->logger->debug("$sequence");
# 					$trigger = 1;
# 				}
# 			}
# 		}
# 		return $self->convertMatchedStringToCoords( $sequence, $seqNameCounter, $start );
# 	}
# 	else {
# 		confess "wrong number of arguments sent to getTrimmedCoords\n";
# 	}
# }

# sub convertMatchedStringToCoords {
# 	my ($self) = shift;

# 	if ( scalar(@_) == 3 ) {
# 		my $sequence              = shift;
# 		my $numberOfAllSequences  = shift;
# 		my $absoluteStartPosition = shift;

# 		my $coordsToReturn;
# 		my $currStart = 0;
# 		my $currEnd   = 0;
# 		my $prevChar  = 0;
# 		my $char      = 0;

# 		for ( my $i = 1 ; $i < length($sequence) ; $i++ ) {
# 			$char = substr( $sequence, $i, 1 );
# 			my $toAdd = 0;

# 			if ( ( $prevChar != $numberOfAllSequences ) && ( $char == $numberOfAllSequences ) ) {
# 				$currStart = $i + $absoluteStartPosition - 1;
# 			}
# 			elsif ( ( $prevChar == $numberOfAllSequences ) && ( $char != $numberOfAllSequences ) ) {
# 				$currEnd = ( $i - 1 ) + ( $absoluteStartPosition - 1 );
# 				$toAdd = 1;
# 			}
# 			elsif ( ( ( $i + 1 ) == length($sequence) ) && ( $char == $numberOfAllSequences ) ) {
# 				$currEnd = $i + $absoluteStartPosition - 1;
# 				$toAdd   = 1;
# 			}

# 			if ($toAdd) {
# 				my $coordToAdd = ',' . $currStart . '..' . $currEnd;
# 				$coordsToReturn .= $coordToAdd;
# 			}
# 			$prevChar = $char;
# 		}
# 		return $coordsToReturn;
# 	}
# 	else {
# 		confess "incorrect number of arguments sent to convertMatchedStringToCoords\n";
# 	}
# }

# sub getSeqNameCoordList {
# 	my ($self) = shift;

# 	if (@_) {
# 		my $hashRef = shift;

# 		my %uniqueNames;
# 		foreach my $query ( keys %{$hashRef} ) {
# 			my $qSeqName = FileInteraction::Fasta::SequenceName->new($query);
# 			$uniqueNames{ $qSeqName->name }->{$query} = $hashRef->{$query};
# 		}
# 		return \%uniqueNames;
# 	}
# 	else {
# 		confess "nothing sent to getSeqNameCoordList\n";
# 	}
# }

# sub isANewMatch {
# 	my ($self) = shift;

# 	if ( scalar(@_) == 3 ) {
# 		my $match         = shift;
# 		my $longestMatch  = shift;
# 		my $longestLength = shift;

# 		#send back as new match if longestLength is at 0
# 		return 1 if ( $longestLength == 0 );

# 		#send back 0 if no match at all
# 		return 0 unless ( length($match) > 0 );

# 		my $longStart;
# 		my $longEnd;
# 		my $matchStart;
# 		my $matchEnd;

# 		if ( $match =~ /\,(\d+)\.\.(\d+)/ ) {
# 			$matchStart = $1;
# 			$matchEnd   = $2;
# 		}
# 		else {
# 			confess "match not correct: $match\n";
# 		}

# 		if ( $longestMatch =~ /\,(\d+)\.\.(\d+)/ ) {
# 			$longStart = $1;
# 			$longEnd   = $2;
# 		}
# 		else {
# 			confess "longestMatch not correct: $longestMatch\n";
# 		}

# 		#we want the best match for each sequenceName
# 		if ( ( $matchStart < $longStart ) || ( $matchEnd > $longEnd ) ) {
# 			return 1;
# 		}
# 		else {
# 			return 0;
# 		}

# 	}
# 	else {
# 		confess "incorrect number of arguments sent to isANewMatch\n";
# 	}
# }

# sub containsAllQueryNames {
# 	my ($self) = shift;

# 	if (@_) {
# 		my $hashRef     = shift;
# 		my $returnValue = 1;
# 		my %sentHashNames;

# 		#get sequence names from hashRef
# 		foreach my $name ( keys %{$hashRef} ) {
# 			my $seqName = FileInteraction::Fasta::SequenceName->new($name);

# 			$self->logger->debug( "containsAllQueryNames: sentName: " . $seqName->name );

# 			$sentHashNames{ $seqName->name } = 1;
# 		}

# 		foreach my $qSeqName ( keys %{ $self->queryNameObjectHash } ) {
# 			unless ( defined $sentHashNames{$qSeqName} ) {
# 				$returnValue = 0;
# 				last;
# 			}
# 		}
# 		return $returnValue;
# 	}
# 	else {
# 		confess "nothing sent to containsAllQueryNames\n";
# 	}
# }

# sub getNumberOfSeqNamesFromHash {
# 	my ($self) = shift;

# 	if (@_) {
# 		my $hashRef = shift;
# 		my %seqNameHash;

# 		foreach my $key ( keys %{$hashRef} ) {
# 			my $seqName = FileInteraction::Fasta::SequenceName->new($key);
# 			$seqNameHash{ $seqName->name } = 1 unless defined $seqNameHash{ $seqName->name };
# 		}
# 		return scalar( keys %seqNameHash );
# 	}
# 	else {
# 		confess "nothing sent to getNumberOfSeqNamesFromHash\n";
# 	}
# }


sub _joinAdjacentNovelRegionsHash {
	my $self = shift;

	if ( scalar(@_) == 1 ) {
		my $hashRef = shift;

		my %joinedHash;

		foreach my $seq ( keys %{$hashRef} ) {
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
	
	foreach my $key(keys %{$self->_comparisonHash}){
		foreach my $key2 (keys %{$self->_comparisonHash->{$key}}){
			$self->logger->debug($key . ' : ' . $key2 . ' : ' . $self->_comparisonHash->{$key}->{$key2});
		}
	}
}




1;
