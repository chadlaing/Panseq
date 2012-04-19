#!/usr/bin/perl

package NovelRegionFinder;

#includes
use strict;
use warnings;
use diagnostics;
use FindBin::libs;
use Carp;
use Mummer::DeltaBlockFactory;
use FileInteraction::Fasta::SequenceName;
use FileInteraction::FileManipulation;
use FileInteraction::Fasta::SequenceRetriever;
#use Mummer::MummerIO;
#use Mummer::MummerParallel;
use Mummer::MummerGPU;
use Pipeline::PanseqShared;

our @ISA=qw/PanseqShared/;

#object creation
sub new{
	my $class=shift;
    my $self = {};
    bless ($self, $class);
    $self->_novelRegionFinderInitialize(@_);   
    return $self;
}

#variables
sub novelRegionFinderMode{
	my $self=shift;
	$self->{'_NRF_novelRegionFinderMode'}=shift // return $self->{'_NRF_novelRegionFinderMode'};
}

sub _skipGatherFiles{
	my $self=shift;
	$self->{'_NRF_skipGatherFiles'}=shift // return $self->{'_NRF_skipGatherFiles'};
}

sub outputFileName{
	my $self=shift;
	$self->{'_NRF_outputFileName'}=shift // return $self->{'_NRF_outputFileName'};
}

sub comparisonHash{
	my $self=shift;
	$self->{'_NRF_comparisonHash'}=shift // return $self->{'_NRF_comparisonHash'};
}

sub _queryFastaNamesHash{
	my $self=shift;
	$self->{'_NRF_queryFastaNamesHash'}=shift // return $self->{'_NRF_queryFastaNamesHash'};
}

sub _adjacentJoiningSize{
	my $self=shift;
	$self->{'_NRF_adjacentJoiningSize'}=shift // return $self->{'_NRF_adjacentJoiningSize'};
}

sub combinedReferenceFile{
	my $self=shift;
	$self->{'_NRF_combinedReferenceFile'}=shift // return $self->{'_NRF_combinedReferenceFile'};
}


#methods
sub _novelRegionFinderInitialize{
	my $self=shift;
	
	#inheritance
	$self->_panseqSharedInitialize(@_);
	
	#anonymous
	$self->_queryFastaNamesHash({});
}

sub validateNovelSettings{
	my($self)=shift;
	
	if(@_){
		my $settingsHashRef=shift;		
		
		foreach my $setting(keys %{$settingsHashRef}){
			my $value = $settingsHashRef->{$setting};
			
			$self->_novel_nucB($self->isAnInt($value)) if $setting eq 'novel_nucB';
			$self->_novel_nucC($self->isAnInt($value)) if $setting eq 'novel_nucC';
			$self->_novel_nucG($self->isAnInt($value)) if $setting eq 'novel_nucG';
			$self->_novel_nucL($self->isAnInt($value)) if $setting eq 'novel_nucL';
			$self->nucD($self->nucDCheck($value)) if $setting eq 'nucD';
			$self->_skipGatherFiles($self->isAnInt($value)) if $setting eq 'skipGatherFiles';
			$self->novelRegionFinderMode($self->novelRegionFinderModeCheck($value)) if $setting eq 'novelRegionFinderMode';	
			$self->combinedQueryFile($value) if $setting eq 'combinedQueryFile';
			$self->combinedReferenceFile($value) if $setting eq 'combinedReferenceFile';
			$self->_adjacentJoiningSize($self->isAnInt($value)) if $setting eq 'adjacentJoiningSize';		
		}
			
		#requirements
		$self->_adjacentJoiningSize(0) unless defined $self->_adjacentJoiningSize;
		$self->missingParam('referenceDirectory') unless (defined $self->referenceDirectory || defined $self->_skipGatherFiles);
		$self->missingParam('novelRegionFinderMode') unless defined $self->novelRegionFinderMode;		
	}
}

sub runNovelRegionFinder{
	my($self)=shift;
	my $configFile=shift;
	
	$self->logger->info("INFO:\tStarting runNovelRegionFinder with configFile: $configFile");
		
	#initialization
	$self->validateNovelSettings($self->getSettingsFromConfigurationFile($configFile));

	#this allows the NovelRegionFinder to be called without having to recombine files
	#and create directories	
	if($self->_skipGatherFiles){
		my $fileManipulator = FileManipulation->new();
		
		$self->logger->info("INFO:\tSkip gathering files. Using previously created combined query file: " . $self->combinedQueryFile);
		$self->queryNameObjectHash($self->getSequenceNamesAsHashRef($fileManipulator->getFastaHeadersFromFile($self->combinedQueryFile)));
	}
	else{
		$self->createDirectories();
		$self->getQueryNamesAndCombineAllInputFiles();
		$self->logger->info("INFO:\tCreated combined query file");
	}	
	
	my $mummer = MummerGPU->new();
	$mummer->run({
		'queryFile'=>$self->combinedQueryFile,
		'referenceFile'=>$self->combinedReferenceFile,
		'mummerDirectory'=>$self->mummerDirectory,
		'deltaFile'=>$self->_baseDirectory . 'combinedDelta.delta',
		'baseDirectory'=>$self->_baseDirectory,
		'numberOfCores'=>$self->_numberOfCores	
	});
	
	$self->logger->info("INFO:\tNucmer finished. Gathering novel regions.");		
	
	#use delta filter to limit size of match 
	#my $newDeltaName = $mummer->deltaFile . '_filtered';
	#my $tempSystemLine = $self->mummerDirectory . 'delta-filter -l 5000 ' . $mummer->deltaFile . ' > ' . $newDeltaName;
	#$mummer->deltaFile($newDeltaName);
		
	#usage is:
	# my $obj->new();
	# $obj->findNovelRegions(<delta file name>, <search type, one of 'common_to_all','unique' and 'no_duplicates'>, 
	#		<hashRef of query names $hash{$name}=1>, <adjacent joining size, or "0" to not join novel regions>);		
		
	my $novelRegionsHashRef = $self->findNovelRegions({
		'deltaFile'=>$mummer->deltaFile,
	});
			
	$self->logger->info("INFO:\tNovel regions gathered. Extracting.");	

	#open output filehandle
	my $outputFileName = $self->_baseDirectory . 'novelRegions.fasta';
	my $novelOutputFH = IO::File->new('>' . $outputFileName) or die "$!";
		
	#order: <novel hash ref>, minimumNovelRegionSize, <output FH if different from STDOUT>
	my $gr = SequenceRetriever->new($self->combinedQueryFile);
	
	$self->logger->debug("DEBUG:\tExtracting and printing regions from Hash
		\tnovelRegionHashRef: $novelRegionsHashRef
		\tnovelOutputFH: $novelOutputFH
		\tcutoffSize: " . $self->minimumNovelRegionSize	
	);
	
	$gr->extractAndPrintRegionsFromHash({
		'novelRegionHashRef'=>$novelRegionsHashRef,
		'novelOutputFH'=>$novelOutputFH	,
		'cutoffSize'=>$self->minimumNovelRegionSize	
	});
	$novelOutputFH->close();
	
	$self->logger->info("INFO:\tNovel regions extracted.");
	return $outputFileName;
}

sub findNovelRegions{
	my($self)=shift;
	
	my $paramsRef=shift;	
	my $deltaFile=$paramsRef->{'deltaFile'} // confess 'deltaFile required in findNovelRegions';		
	
	#construct hash of all queryNameObjectHash fasta headers
	$self->_queryFastaNamesHash($self->_getAllFastaHeadersFromNameHash());
	
	#build hash of all comparisons
	$self->buildHashOfAllComparisons($deltaFile);
		
	#get the novel regions
	return $self->getNovelRegions();	
}

sub _getAllFastaHeadersFromNameHash{
	my $self=shift; 
	
	#get hash of query fasta headers for novelRegionFinder
	my %fastaNamesHash;
	my $sr = SequenceRetriever->new($self->combinedQueryFile);
	
	foreach my $seqName(keys %{$self->queryNameObjectHash}){
		my $headersRef = $self->queryNameObjectHash->{$seqName}->arrayOfHeaders;
		foreach my $header(@{$headersRef}){
			my $headerLength = length($sr->extractRegion($header));			
			$fastaNamesHash{$header}=$headerLength;
			
			$self->logger->debug("DEBUG:\tPopulating fasta header hash $header:$headerLength");
		}
	}
	return \%fastaNamesHash;
}

sub getNovelRegions{
	my($self)=shift;
	
	if(defined $self->novelRegionFinderMode){
	
		my $novelCoordsHashRef;
		
		if($self->novelRegionFinderMode eq 'no_duplicates'){
			$novelCoordsHashRef= $self->getNoDuplicates();	
		}
		elsif($self->novelRegionFinderMode eq 'common_to_all'){
			$novelCoordsHashRef = $self->getCommonToAll();
		}
		elsif($self->novelRegionFinderMode eq 'unique'){
			$novelCoordsHashRef = $self->getUnique();
		}
		else{
			confess "incorrect type sent to getNovelRegions\n";
		}
		return $self->joinAdjacentNovelRegionsHash($novelCoordsHashRef);		
	}
	else{
		confess "comparisonType not defined in getNovelRegions\n";
	}
}

sub getUnique{
	my($self)=shift;	

	my %uniqueMatches;
		
	foreach my $query(keys %{$self->comparisonHash}){
		foreach my $ref(keys %{$self->comparisonHash->{$query}}){
			
			$self->logger->debug('DEBUG:' . "\tGet unique: Query: $query Ref: $ref");
			
			next if $query eq $ref;
			$uniqueMatches{$query} .= $self->comparisonHash->{$query}->{$ref};
		}
	}
	my $sortedHashRef = $self->sortAndCompileCoordsHash(\%uniqueMatches);
	return $self->getNegativeImageCoords($sortedHashRef);
}

sub getCommonToAll{
	my($self)=shift;
	
	my %commonCoordsHash;
		
	#use one query sequence as the one from which to extract the sequence data
	#common in all allows this
	#for each query, check that all other queries have a hit
	#using "one" as the w.r.t sequence, get all reference hits for the query
	#take negative image
				
	my %oneSeqHash; #this is a hash so that it can be fed into sortAndCompileCoordsHash
	
	my @queryNames = keys %{$self->queryNameObjectHash};
	my $oneSeq = SequenceName->new($queryNames[0]);	
	
	$self->logger->debug("DEBUG:\tCTA: oneSeq: ". $oneSeq->name);
	
	foreach my $query(keys %{$self->comparisonHash}){
		$self->logger->debug("DEBUG:\tCTA: query: $query");
		
		my $querySeqName = SequenceName->new($query);
		next unless $querySeqName->name eq $oneSeq->name;
		
		#get the "oneSeq" coords
		foreach my $hit(keys %{$self->comparisonHash->{$query}}){	
			$self->logger->debug("DEBUG:\tCTA: hit: $hit");			
			
			if (defined $self->_queryFastaNamesHash->{$hit}){
				next;
				#we only want sequence matches in the reference, for the negative image				
			}	
			else{
				$self->logger->debug("DEBUG:\t$hit not defined in queryFastaNamesHash");
			}	
			$oneSeqHash{$query} .= $self->comparisonHash->{$query}->{$hit};
			
			$self->logger->debug("DEBUG:\tCTA: oneSeq: " . $oneSeq->name . " oneSeqHash{query}:" . $oneSeqHash{$query});		
		}	
	}

	return $self->trimOneSeqToAllQuerySequences(
		$self->getNegativeImageCoords(
			$self->sortAndCompileCoordsHash(\%oneSeqHash)
		)	
	);	
	
}

sub trimOneSeqToAllQuerySequences{
	my($self)=shift;
	
	if(@_){
		my $hashRef=shift; #this is the hashRef of all reference regions not found in oneSeq
		my %trimmedHash;
		
		$self->logger->debug("DEBUG:\tTrim: begin trim, number of hashRef keys is: " . scalar keys %{$hashRef});
		
		#should only be one sequence name, but possibly multiple keys due to the oneSeq selection in getCommonToAll
		foreach my $oneSeqQueryName(keys %{$hashRef}){		
			my $comparisonHashRef = $self->comparisonHash->{$oneSeqQueryName};
			
			#the following loop stores only query hits in the %coords hash
			my %coords;
			foreach my $comparisonHit(keys %{$comparisonHashRef}){
				$self->logger->debug("DEBUG:\tTrim: comparisonHit: $comparisonHit");	
				next unless defined $self->_queryFastaNamesHash->{$comparisonHit};				
				$self->logger->debug("DEBUG:\tTrim: made the cut: $comparisonHit: ". $comparisonHashRef->{$comparisonHit});				
				$coords{$comparisonHit} = $comparisonHashRef->{$comparisonHit};				
			}			
			$trimmedHash{$oneSeqQueryName} = $self->getCommonCoordsWRTinputCoords($hashRef->{$oneSeqQueryName},\%coords); #input: coords as string , hashRefToCheck
		}
		return \%trimmedHash;
	}
	else{
		confess "nothing sent to trimOneSeqToAllQuerySequences\n";
	}
}

sub getCommonCoordsWRTinputCoords{
	my($self)=shift;
	
	if(scalar(@_)==2){
		my $inputCoordString = shift // die "undefined inputCoordString in getCommonCoordsWRTinputCoords\n";
		my $hashRef = shift; #all the query seqs and their hits (no ref hits)
		my $newString;
		
		$self->logger->debug("DEBUG:\tWRT inputCoordString: $inputCoordString");
		
		#algorithm
		#get seqName list for each inputCoord (all query sequence matches for the given oneSeq are in the hashRef)
		#cycle through seqName list, trimming as it goes
		
		my $seqNameCoordsHashRef = $self->getSeqNameCoordList($hashRef);
		
		my $iteration=0;

		while($inputCoordString =~ /(\,\d+\.\.\d+)/gc){
			my $currentCoords = $1;
			
			$self->logger->debug("DEBUG:\tWRT currentCoords: $currentCoords");
			
			my $returnedValue = $self->getTrimmedCoords($currentCoords, $seqNameCoordsHashRef);
			
			if (defined $returnedValue){
				 $self->logger->debug("DEBUG:\tWRT afterTrim: $returnedValue");
			}
			else{
				$self->logger->debug("DEBUG:\tWRT afterTrim: not common to all");
				next;
			}
			
			$newString .=$returnedValue;			
			$self->logger->debug("DEBUG:\tWRT newString: $newString");						
		}
		return $newString;		
	}
	else{
		confess "nothing sent to getCommonCoordsWRTinputCoords\n";
	}
}

sub getTrimmedCoords{
	my($self)=shift;
	
	#This sub creates a string of 0s representing the length of coords sent to the sub 0000000000
	#It then increments each position by 1 for every match eg. 0001111100, then 00012222200
	#and outputs only the region from the original string that has numbers equal the total number of query sequences	
	
	
	
	if(scalar(@_)==2){
		my $coords=shift;
		my $hashRef=shift; # $hashRef->{seqName}->{seqFastaHeader}=<coords>
		
		#get start/stop coords
		my $start;
		my $end;
		
		if($coords =~ /\,(\d+)\.\.(\d+)/){
			$start = $1;
			$end = $2;
		}
		else{
			confess "no coords sent to getTrimmedCoords\n";
		}
		
		my $offset = $start-1;
		
		$self->logger->debug("DEBUG:\tGTC coords: $coords start: $start end: $end offset:$offset");
		
		my $sequence = '0' . ('0' x ($end - $start +1)); #to account for 1 offset
		
		#check that each seqName has a hit, and where it is
		my $seqNameCounter=0;
		my $trigger=1;
		my $previousName;
		
		foreach my $sName(keys %{$hashRef}){
			if($trigger==1){
				$seqNameCounter++;
				$self->logger->info("INFO:\tMatching $sName");
			}
			else{
				$self->logger->warn("WARN:\tCommon sequence not present in getTrimmedCoords for $previousName");
				last;
			}			
			
			$trigger=0;
			$previousName=$sName;
			
			$self->logger->debug("DEBUG:\tsName: $sName seqNameCounter: $seqNameCounter");
			
			foreach my $eachMatch(keys %{$hashRef->{$sName}}){
				my $eachCoord = $hashRef->{$sName}->{$eachMatch};
				
				$self->logger->debug("DEBUG:\teachMatch $eachMatch");
				
				while($eachCoord =~ /\,(\d+)\.\.(\d+)/gc){
					
					my $eachStart=$1;
					my $eachEnd=$2;						
					
					$self->logger->debug("DEBUG:\t$eachMatch $eachStart:$eachEnd");
					
					#check to make sure each start/end is within the range
					if($eachStart > $end){
						$self->logger->debug("DEBUG:\tNEXT eachStart greater than end");
						next;
					}
					
					if($eachEnd < $start){
						$self->logger->debug("DEBUG:\tNEXT eachEnd less than start");
						next;						
					}
					$self->logger->debug("DEBUG:\tMATCH eachStart: $eachStart eachEnd: $eachEnd");
					$eachStart = $start if $eachStart < $start; 
					$eachEnd = $end if $eachEnd > $end;
										
					my $eachLength = $eachEnd - $eachStart +1;
					
					$self->logger->debug("DEBUG:\tBefore adjustment eachStart: $eachStart start:$start eachEnd: $eachEnd end:$end");
					
					#account for difference between absolute and relative locations
					$eachStart = $eachStart - $offset;
					$eachEnd = $eachEnd - $offset;
					
					#get substring, increase count by 1, replace in sequence string
					$self->logger->debug("DEBUG:\tAfter adjustment eachStart: $eachStart eachEnd: $eachEnd");
					
					my $substring = substr($sequence,$eachStart, ($eachEnd - $eachStart +1));
					my $prevNum = $seqNameCounter -1;

					$substring =~ s/$prevNum/$seqNameCounter/g;
					substr($sequence,$eachStart,$eachLength)=$substring;
					$self->logger->debug("$sequence");
					$trigger=1;
				}			
			}
		}
		return $self->convertMatchedStringToCoords($sequence,$seqNameCounter, $start);		
	}
	else{
		confess "wrong number of arguments sent to getTrimmedCoords\n";
	}
}

sub convertMatchedStringToCoords{
	my($self)=shift;
	
	if(scalar(@_)==3){
		my $sequence =shift;
		my $numberOfAllSequences=shift;
		my $absoluteStartPosition=shift;
		
		my $coordsToReturn;
		my $currStart=0;
		my $currEnd=0;
		my $prevChar=0;
		my $char=0;

		for(my $i=1; $i < length($sequence); $i++){
			$char = substr($sequence,$i,1);
			my $toAdd=0;

			if(($prevChar != $numberOfAllSequences) && ($char == $numberOfAllSequences)){
				$currStart=$i + $absoluteStartPosition -1;
			}
			elsif(($prevChar == $numberOfAllSequences) && ($char != $numberOfAllSequences)){
				$currEnd = ($i-1) + ($absoluteStartPosition-1);
				$toAdd=1;
			} 
			elsif((($i+1) == length($sequence)) && ($char == $numberOfAllSequences))
			{
				$currEnd = $i + $absoluteStartPosition -1;
				$toAdd=1;
			}
			
			if($toAdd){
				my $coordToAdd = ',' . $currStart . '..' . $currEnd;
				$coordsToReturn .= $coordToAdd;
			}
			$prevChar=$char;
		}
		return $coordsToReturn;
	}
	else{
		confess "incorrect number of arguments sent to convertMatchedStringToCoords\n";
	}
}


sub getSeqNameCoordList{
	my($self)=shift;
	
	if(@_){
		my $hashRef=shift;
		
		my %uniqueNames;
		foreach my $query(keys %{$hashRef}){
			my $qSeqName = SequenceName->new($query);
			$uniqueNames{$qSeqName->name}->{$query}=$hashRef->{$query};
		}
		return \%uniqueNames;
	}
	else{
		confess "nothing sent to getSeqNameCoordList\n";
	}
}

sub isANewMatch{
	my($self)=shift;
	
	if(scalar(@_)==3){
		my $match=shift;
		my $longestMatch=shift;
		my $longestLength=shift;
		
		#send back as new match if longestLength is at 0
		return 1 if($longestLength ==0);
		
		#send back 0 if no match at all
		return 0 unless(length($match) > 0);
		
		my $longStart;
		my $longEnd;
		my $matchStart;
		my $matchEnd;
		
		if($match =~ /\,(\d+)\.\.(\d+)/){
			$matchStart=$1;
			$matchEnd=$2;
		}
		else{
			confess "match not correct: $match\n";
		}
		
		if($longestMatch =~ /\,(\d+)\.\.(\d+)/){
			$longStart=$1;
			$longEnd=$2;
		}
		else{
			confess "longestMatch not correct: $longestMatch\n";
		}
		
		#we want the best match for each sequenceName
		if(($matchStart < $longStart) || ($matchEnd > $longEnd)){
			return 1;
		}
		else{
			return 0;
		}
		
	}
	else{
		confess "incorrect number of arguments sent to isANewMatch\n";
	}
}


sub containsAllQueryNames{
	my($self)=shift;
	
	if(@_){
		my $hashRef=shift;
		my $returnValue=1;
		my %sentHashNames;
		
		#get sequence names from hashRef
		foreach my $name(keys %{$hashRef}){
			my $seqName = SequenceName->new($name);
			
			$self->logger->debug("DEBUG:\tcontainsAllQueryNames: sentName: " . $seqName->name);
			
			$sentHashNames{$seqName->name}=1;
		}
		
		foreach my $qSeqName(keys %{$self->queryNameObjectHash}){
			unless(defined $sentHashNames{$qSeqName}){
				$returnValue=0;
				last;
			}
		}		
		return $returnValue;
	}
	else{
		confess "nothing sent to containsAllQueryNames\n";
	}
}

sub getNumberOfSeqNamesFromHash{
	my($self)=shift;
	
	if(@_){
		my $hashRef=shift;
		my %seqNameHash;
		
		foreach my $key(keys %{$hashRef}){
			my $seqName = SequenceName->new($key);
			$seqNameHash{$seqName->name}=1 unless defined $seqNameHash{$seqName->name};
		}
		return scalar(keys %seqNameHash);
	}
	else{
		confess "nothing sent to getNumberOfSeqNamesFromHash\n";
	}
}

sub getNoDuplicates{
	my($self)=shift;
	

	my %nonDuplicatedCoords;
		
	$self->logger->info("INFO:\tGetting no duplicates");
		
	#algorithm
	#A vs refs
	#B vs A + refs
	#C vs. A + B + refs
	#...	
		
	my %checkedQuery;
	
	foreach my $query(keys %{$self->comparisonHash}){
		next unless defined $self->_queryFastaNamesHash->{$query};
		
		$self->logger->debug("DEBUG:\tGET_NO_DUPLICATES: query: $query");
		foreach my $ref(keys %{$self->comparisonHash->{$query}}){
			next if $query eq $ref;
				
			#if a query name as a reference hit
			if(defined $self->_queryFastaNamesHash->{$ref}){
				next unless ((defined $checkedQuery{$ref}));
			}	
					
			$nonDuplicatedCoords{$query} .= $self->comparisonHash->{$query}->{$ref};
			$self->logger->debug("DEBUG:\tGET_NO_DUPLICATES: ref: $ref non_dup_coords: $nonDuplicatedCoords{$query}");	
		}
		$checkedQuery{$query}=1;
	}
	$self->logger->info("INFO:\tSorting and compiling novel regions");
	my $sortedHashRef = $self->sortAndCompileCoordsHash(\%nonDuplicatedCoords);
		
	$self->logger->info("INFO:\tGetting negative image of coords");
	return $self->getNegativeImageCoords($sortedHashRef);		
}

sub getNegativeImageCoords{
	my($self)=shift;
	
	if(@_){
		my $hashRef=shift;
		my %negativeHash;
		#hash of $hash{$query}=<coords in order>
		
		foreach my $key(keys %{$hashRef}){
			my $missingCoords;
			my $prevEnd=0;
			
			#get the negative image coords
			my $nextKey = $hashRef->{$key};
			while($nextKey =~ /\,(\d+)\.\.(\d+)/gc){
				my $start=$1;
				my $end=$2;
				
				$self->logger->debug("DEBUG:\tNEG:\tMatch coords: start:$start end:$end");
				
				if(($start > ($prevEnd+1)) && ($start != 1) && (($start +1) < $end)){
					
					#kludge for 1..1 init match
					my $advancement=1;
					if($prevEnd eq '1'){
						$advancement=0;
					}	
					$missingCoords .= ',' . ($prevEnd + $advancement) . '..' . ($start -1);
					
					$self->logger->debug("DEBUG:\tNEG:\tpassed checks: missingCoords: $missingCoords prevEnd:$prevEnd start:$start");
				}				
				$prevEnd=$end;				
			}
			
			#check the last end
			$self->logger->debug("DEBUG:\tKey: $key");
			if(($self->_queryFastaNamesHash->{$key} > $prevEnd) && ($prevEnd != 0)){
				
				$self->logger->debug("DEBUG:\tNEG: prevEnd: $prevEnd");
				
				#kludge for 1..1 init match
				my $advancement=1;
				if($prevEnd eq '1'){
					$advancement=0;
				}				
				$missingCoords .= ',' . ($prevEnd +$advancement) . '..' . ($self->_queryFastaNamesHash->{$key});
			}
						
			$self->logger->debug("DEBUG:\tNEG: original of $key: $nextKey");
			if(defined $missingCoords){
				$self->logger->debug("DEBUG:\tNEG: negative: $missingCoords");
			}
			else{
				$self->logger->debug("DEBUG:\tNEG: no missing coords");
			}		
			
			$negativeHash{$key}=$missingCoords if defined $missingCoords;
		}
		return \%negativeHash;
	}
	else{
		confess "nothing sent to getNegativeImageCoords\n";
	}
}


sub sortAndCompileCoordsHash{
	my($self)=shift;
	
	if(@_){
		my $hashRef=shift;
		
		foreach my $query(keys %{$hashRef}){
			
			$self->logger->debug('DEBUG:' . "\tSorting and compiling $query");
			
			$hashRef->{$query} = $self->sortAndCompileCoordString($hashRef->{$query});
		}
		return $hashRef;
	}
	else{
		confess "nothing sent to sortAndCompileCoordsHash\n";
	}
	
}

sub sortAndCompileCoordString{
	my($self)=shift;
	
	if(@_){
		my $coords=shift;
		my %coordsHash;
		my $newString;
		
		#create hash with first coord as key
		
		
		while($coords =~ /(\,(\d+)\.\.\d+)/gc){
			#do the sorting and compiling
			my $coord=$1;
			my $startKey = $2;
			
			#account for duplicate start coords --> otherwise problems
			while(defined $coordsHash{$startKey}){
				$startKey++;
			}			
			$coordsHash{$startKey}=$coord;
		}
		
		#sort hash, reconstruct string
		foreach my $key(sort {$a <=> $b} keys %coordsHash){
			if(defined $newString){
				my $start;
				my $end;
				my $sortCoord;
				if($coordsHash{$key} =~ /(\,(\d+)\.\.(\d+))/){
					$sortCoord=$1;
					$start=$2;
					$end=$3;
					
					$self->logger->debug('DEBUG:' . "\tSort and compile: Before update: Key: $key start: $start end: $end");
					
				}
				else{
					confess $coordsHash{$key} . ' did not match!' . "\n";
				}
				$newString = $self->updateAlignmentCoords($newString,$start,$end);
				$self->logger->debug('DEBUG:' . "\tSort and compile: After update: Key: $key start: $start end: $end");

			}
			else{
				$newString = $coordsHash{$key};
				
			}		
		}
		$self->logger->debug('DEBUG:' . "\tSort and compile: New string: $newString");
		return $newString;
	}
	else{
		confess "nothing sent to sortAndCompileCoordString\n";
	}
}

sub joinAdjacentNovelRegionsHash{
	my($self)=shift;
	
	if(scalar(@_)==1){
		my $hashRef=shift;
		
		my %joinedHash;
		
		foreach my $seq(keys %{$hashRef}){ 
			unless(defined $hashRef->{$seq}){
				$self->logger->debug("DEBUG:\t$seq has no novel coords");
				next;
			}
			$joinedHash{$seq}= $self->joinAdjacentNovelRegionsString($hashRef->{$seq});
		}
		return \%joinedHash;
	}
	else{
		confess ("incorrect number of arguments sent to joinAdjacentNovelRegions\n");
	}
}

sub joinAdjacentNovelRegionsString{
	my($self)=shift;
	
	if(scalar(@_)==1){
		my $coordsString=shift // confess ("UNDEF coordsString sent to joinAdjacentNovelRegionsString");
		
		$self->logger->debug("DEBUG:\tcoordsString in joinAdjacentNovelRegionsString is: $coordsString");		
		return $coordsString if ($self->_adjacentJoiningSize == '0');
		
		my $prevStart;
		my $prevEnd;
		my $start;
		my $end;
		my $prevCoords;
		my $currCoords;
		my $newCoordsString = $coordsString;
		
		my $count=0;
		while($coordsString =~ /(\,(\d+)\.\.(\d+))/gc){

			$currCoords = $1;
			$start=$2;
			$end=$3;
						
			next unless defined $prevStart;
			
			if(($start - $prevEnd) < $self->_adjacentJoiningSize){
				my $oldCoords = $prevCoords . $currCoords;
				my $newCoords = ',' . $prevStart . '..' . $end;
				$newCoordsString =~ s/$oldCoords/$newCoords/;
			}			
		}
		continue{
			$prevStart=$start;
			$prevEnd=$end;
			$prevCoords=$currCoords;
		}
		return $newCoordsString;
	}
	else{
		confess "incorrect number of arguments sent to joinAdjacentNovelRegionsString\n";
	}
}

sub buildHashOfAllComparisons{
	my($self)=shift;
	
	$self->logger->info("INFO:\tBuilding hash of all comparisons");
	
	#create a self vs self entry to accomodate regions that are completely novel and have no hit
	#unless($self->novelRegionFinderMode eq 'common_to_all'){
		foreach my $fastaHeader(keys %{$self->_queryFastaNamesHash}){
			$self->addToComparisonHash($fastaHeader,'_initialize',',1..1');
		}
	#}
	
	if(@_){
		my $inputFileName = shift;
		
		my $inputFH = IO::File->new('<' . $inputFileName) or die "$!";
		my $dbFactory = DeltaBlockFactory->new($inputFH);
		
		while(my $block = $dbFactory->nextDeltaBlock(1)){	#turns on absolute start/end positions with the true value
			$self->updateComparisonHash($block,'query');
			$self->updateComparisonHash($block,'reference'); #this allows for running a single Nucmer comparison but still getting the all vs. all output
		}
		$inputFH->close();		
	}
	else{
		confess "nothing sent to buildHashOfAllComparisons\n";
	}
}

sub updateComparisonHash{
	my($self)=shift;
	
	if(scalar(@_)==2){
		my $block=shift;
		my $type=shift; #can be either query or reference, and determines how the hash is built
						#'query' means $hashRef->{query}->{ref}=<query coords>
						#'reference' means $hashRef->{ref}->{query}=<ref coords>
						#also determines what length value is sent to comparisonLengthHash
		my $alignmentCoords; 
		
		$self->logger->debug("DEBUG:\tUCH begin");
		
		my $key1;
		my $key2;
		my $length;
		my $startCoord;
		my $endCoord;
		
		if($type eq 'query'){
			$key1=$block->queryName;
			$key2=$block->refName;
			$length=$block->queryLength;
			$startCoord=$block->queryStart;
			$endCoord=$block->queryEnd;
		}
		elsif($type eq 'reference'){
			$key1=$block->refName;
			$key2=$block->queryName;
			$length=$block->refLength;
			$startCoord=$block->refStart;
			$endCoord=$block->refEnd;
		}
		else{
			confess "incorrect comparison type in updateComparisonHash\n";
			exit(1);
		}
		
		if($startCoord > $endCoord){
			$self->logger->fatal("FATAL\tstartCoord:$startCoord bigger than endCoord:$endCoord");
			exit(1);			
		}
		
		if((defined $self->comparisonHash) && (defined $self->comparisonHash->{$key1}->{$key2})){
			$self->logger->debug("DEBUG:\tUCH: Alignment coords previously exist for $key1:$key2");
			$alignmentCoords = $self->comparisonHash->{$key1}->{$key2};
		}
		
		if(defined $alignmentCoords){
			$self->logger->debug("DEBUG:\tUCH:Before update: $key1:$key2 $alignmentCoords");
			$self->addToComparisonHash($key1, $key2, 
				$self->updateAlignmentCoords($alignmentCoords,  $startCoord, $endCoord));
		}
		else{
			$alignmentCoords = ',' . $startCoord . '..' . $endCoord;
			
			$self->logger->debug("DEBUG:\tUCH:New alignment: $key1:$key2 $alignmentCoords");
			
			$self->addToComparisonHash($key1, $key2, $alignmentCoords);
		}	
	}
	else{
		confess "nothing sent to updateComparisonHash or comparisonType undefined\n";
	}
}

sub updateAlignmentCoords{
	my($self)=shift;
	
	if(scalar(@_)==3){
		my $alignmentCoords =shift;
		my $qStart=shift; #current values 'q', checking whether they can modify the start/end of any coords
		my $qEnd=shift;
		
		$self->logger->debug('DEBUG:' . "\tUAC: qStart: $qStart qEnd: $qEnd");
		
		my $newStart;
		my $newEnd;
		my $prevStart;
		my $prevEnd;
		my $lastTime=0;
		my $duplicate=0;
		my $oldString;
		my $newString;

			while($alignmentCoords =~ /\,(\d+)\.\.(\d+)/gc){
				$prevStart = $1;
				$prevEnd =$2;
				
				$self->logger->debug('DEBUG:' . "\tUAC: prevStart: $prevStart prevEnd: $prevEnd");
				
				#check for a q start that is less than or between the current values
				if($qStart > $prevEnd){
					next;
				}
				
				#if both qStart and qEnd are within the current set, no need to add the duplicate
				if(($qStart >= $prevStart) && ($qEnd <= $prevEnd)){
					$duplicate=1;
					last;
				}
				
				#initial values
				$newStart=$prevStart;
				$newEnd=$prevEnd;
				
				#extend the beginning range if needed
				if(($qStart < $prevStart) && ($qEnd >= $prevStart)){
					$newStart = $qStart;
					$lastTime=1;
				}
				
				#extend the end range if needed
				if(($qEnd > $prevEnd) && ($qStart <= $prevEnd)){
					$newEnd = $qEnd;
					$lastTime=1;
				}
				
				if($lastTime){
					$oldString = ',' . $prevStart . '..' . $prevEnd;
					last;
				}			
			}#end while
		
		
		#double check something has changed, then update
			if(defined $oldString){	
				$self->logger->debug('DEBUG:' . "\tUAC: Old string found: newStart: $newStart newEnd: $newEnd");
				
				$newString = ',' . $newStart . '..' . $newEnd;		
				$alignmentCoords =~ s/\Q$oldString\E/$newString/;
			}
			elsif(!$duplicate){
				#this means nothing was modified, but the new values need to be added
				
				$self->logger->debug('DEBUG:' . "\tUAC: New string found: qStart: $qStart qEnd: $qEnd");
				
				$newString = ',' . $qStart . '..' . $qEnd;
				$alignmentCoords .= $newString;
			}
			#otherwise, duplicate value, nothing added	
		$self->logger->debug('DEBUG:' . "\tUAC: alignmentCoords: $alignmentCoords");
		return $alignmentCoords;		
	}
	else{
		confess "incorrect number of arguments sent to updateAlignmentCoords\n";
	}
}

sub addToComparisonHash{
	my($self)=shift;
	
	if(scalar(@_)==3){
		my $q=shift;
		my $r=shift;
		my $coords=shift;
		
		if(defined $self->comparisonHash){
			$self->comparisonHash->{$q}->{$r}=$coords;
		}
		else{
			my %tempHash;
			$tempHash{$q}->{$r}=$coords;
			$self->comparisonHash(\%tempHash);
		}
	}
	else{
		confess "wrong number of arguments sent to addToComparisonHash\n";
	}
}


sub novelRegionFinderModeCheck{
	my($self)=shift;
	
	if(@_){
		my $mode=shift;
		
		if(($mode eq 'no_duplicates') || ($mode eq 'common_to_all') || ($mode eq 'unique')){
			return $mode;
		}
		else{
			confess "$mode is not a valid novelRegionFinderMode!\n",
				"valid modes are no_duplicates, common_to_all and unique.\n";
		}
	}
	else{
		confess "nothing sent to novelRegionFinderModeCheck\n";
	}
}


1;