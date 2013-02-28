#!/usr/bin/env perl
package Modules::Alignment::SNPFinder;

use FindBin;
use lib "$FindBin::Bin/../../";
use Modules::Fasta::SequenceName;
use Log::Log4perl;

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_initialize(@_);
    return $self;
}

#class vars
sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}

sub allowableChars{
	my $self=shift;
	$self->{'_allowableChars'} = shift // return $self->{'_allowableChars'};
}

sub orderedNames{
	my $self=shift;
	$self->{'_queryNameOrderHash'}=shift // return $self->{'_queryNameOrderHash'};
}


sub alignedFastaSequences{
	my $self=shift;
	$self->{'_alignedFastaSequences'}=shift // return $self->{'_alignedFastaSequences'};
}

sub resultNumber{
	my $self=shift;
	$self->{'_resultNumber'}=shift // return $self->{'_resultNumber'};
}

sub startBpHashRef{
	my $self=shift;
	$self->{'_startBpHashRef'}=shift // return $self->{'_startBpHashRef'};
}

#methods
sub _initialize{
	my($self)=shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::SNPFinder\n");
	
	my %params = @_;
    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::Alignment::SNPFinder");
		}
	}
	
	#defaults		
	$self->allowableChars({
		'A'=>1,
		'C'=>1,
		'T'=>1,
		'G'=>1,
		'a'=>1,
		'c'=>1,
		't'=>1,
		'g'=>1
	});
}

sub findSNPs{
	my($self)=shift;
	
	my ($alignmentLength,$alignedHashRef) = $self->_getHashOfFastaAlignment();

	my @orderedResults; 
	my $addedNames=0;
	for my $position(0..($alignmentLength-1)){
		my $resultLine = $self->_getSingleBaseResult($position,$alignedHashRef);

		#add the contig names as tab-delimited values at the start of each new contig
		if(defined $resultLine && $addedNames==0){
			$addedNames=1;
			$resultLine .= $self->_getFastaNamesForOutput($alignedHashRef);
		}

		if(defined $resultLine){
			push @orderedResults,$resultLine;
		}
	}
	
	if(scalar(@orderedResults) > 0){
		return \@orderedResults;
	}
	else{
		return undef;
	}
	
}

=head3 _getFastaNamesForOutput

Returns an ordered list of contig names present in the current result for each genome.
The final output is as below:
name:data:position:contigNames
eg.
	
snp_1000000012	A 			A 			T 			A  			100004 		100004 			544 		 	100005	Acontig001 	Bcontig01004 	Ccontig0034 	Dcontig000043

=cut

sub _getFastaNamesForOutput{
	my $self=shift;
	my $alignedHashRef=shift;

	my $fastaNameLine='';
	foreach my $sName(@{$self->orderedNames}){
		#it could be there were no BLAST results for a given name
		#in that case, add 'N/A'
		if(defined $alignedHashRef->{$sName}->{'fasta'}){
			$fastaNameLine .=("\t" . $alignedHashRef->{$sName}->{'fasta'});
		}
		else{
			$self->logger->debug("$sName not found in alignedHashRef");
			$fastaNameLine .=("\t" . 'N/A');
		}
	}
	return $fastaNameLine;
}


=head3

Given a hash of the fasta sequence alignment,
and the position of the alignment, checks if a SNP exists.
SNPs must be of $self->allowableChars.
Returns a tab-delimited line of results, in the order of
$self->orderedNames.

=cut


sub _getSingleBaseResult{
	my $self = shift;
	my $position=shift;
	my $alignedHashRef=shift;

	my %baseTypes;
	my $baseLine='';
	my $positionLine='';

	foreach my $name(@{$self->orderedNames}){
		my $base;
		#in the case where there was no BLAST hit, there is no hash key to look up
		#need to fill in the position as '-'
		if(exists $alignedHashRef->{$name}->{'fasta'}){
			$base = substr($alignedHashRef->{$name}->{'sequence'},$position,1); 

			unless(defined $base){
				$self->logger->warn("name: $name\nfasta: " . $alignedHashRef->{$name}->{'fasta'}. "\nseq: " . $alignedHashRef->{$name}->{'sequence'} . "\npos: $position");
			}

			if(defined $self->allowableChars->{$base}){
				$baseTypes{$base}=1;
			}
			
			$baseLine .= ("\t$base");	
			my $startBp = $self->startBpHashRef->{$alignedHashRef->{$name}->{'fasta'}};

			$positionLine .= ("\t" . ($startBp + $position));			
		}
		else{
			$baseLine .= ("\t" . '-');
			$positionLine .= ("\t" . '-');
		}
	}

	if(scalar keys %baseTypes > 1){
		return (('snp_' . ($self->resultNumber + $position)) . $baseLine . $positionLine);
	}
	else{
		return undef;
	}
}


=head3

Given the FASTA alignment produced by Muscle, create a hash where the
name is based on the Modules::Fasta::SequenceName->name and
$hashRef->{'name'}->{'fasta'}="fasta header"
$hashRef->{'name'}->{'sequence'}="DNA sequence"
the {'fasta'} key contains the fasta header,
and the {'sequence'} key contains the DNA sequence.
Also compute the alignment length and return it as the first argument of the two argument return.

=cut

sub _getHashOfFastaAlignment{
	my $self = shift;

	my $results={};
	my $name;
	my $alignmentLength;

	my $numberOfLines = scalar @{$self->alignedFastaSequences};

	foreach my $line(@{$self->alignedFastaSequences}){
		$line =~ s/\R//g;
		
		if($line =~ /^>(.+)/){		
			my $sn = Modules::Fasta::SequenceName->new($line);
			$name = $sn->name;

			#we don't need or want the '>'; all names are stored without the fasta header signifier
			$results->{$name}->{'fasta'}=$1;
		}
		else{
			if(defined $results->{$name}->{'sequence'}){
				$results->{$name}->{'sequence'}=$results->{$name}->{'sequence'} . $line;	
			}
			else{
				$results->{$name}->{'sequence'}=$line;
			}
					
		}
	}
	$alignmentLength=length($results->{$name}->{'sequence'});
	return ($alignmentLength,$results);
}
1;
