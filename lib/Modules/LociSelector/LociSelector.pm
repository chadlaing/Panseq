#!/usr/bin/env perl

=pod

=head1 Modules::LociSelector::LociSelector

Modules::LociSelector::LociSelector - Heuristically finds loci that offer maximum discrimination among data.

=head1 SYNOPSIS


    use Modules::LociSelector::LociSelector;

    my $finder = Modules::LociSelector::LociSelector->new(
        inpurFile=>'filename',
        lociNumver=>'best' | number,
        maximizePod => 0  #optional, default is 0 
    );
    $finder->run();

=head1 DESCRIPTION

Modules::LociSelector::LociSelector constructs loci sets that are maximized with respect to the unique number of fingerprints produced among the input sequences as well as the discriminatory power of the loci among the input sequences. The final loci set is iteratively built, in the following steps, given a tab-delimited table with loci names in the first column, sequence names in the first row, and single character data filling the matrix. Missing data is denoted by the characters '?', '-', or '.' :

(1) Each potential available locus is evaluated for the number of unique fingerprints that would result from its addition to the final loci set. All loci that would generate the maximum number of unique fingerprints in this respect are evaluated in step (2).

(2) All loci from step (1) are evaluated for their discriminatory power among the sequences, which is given as points of discrimination (POD). The POD for a locus is calculated as follows.

A listing of all possible pair-wise comparisons is constructed; for example, if the input table consisted of three sequences, A, B and C, the list would consist of A-B, A-C and B-C. Next, it is determined whether or not the sequences in each pair-wise comparison contain the same single character denoting the locus state. If they do, a value of 0 is assigned; if they differ a value of 1 is assigned. The POD is then the summation of all pair-wise comparisons that differ for that locus. With our previous example, if A-B = 1, A-C = 1 and B-C = 0, the POD for that locus would be 2.

(3) The locus with the highest value from step (2) is selected for addition to the final loci set and removed from the pool of candidate loci. If two or more loci tie in value, one is randomly selected. If all possible unique fingerprints have been found, the algorithm continues with (4); if additional unique fingerprints are possible, the algorithm continues with (5).

(4) Sequence pairs for which the allele of the locus chosen in (3) differ are temporarily excluded from the analysis ("masked"). This ensures loci that differ between other pairs of strains are preferentially considered. Consider our A, B and C example with pair-wise comparisons of A-B = 1, A-C = 1 and B-C = 0. In the case of this locus being chosen, the sequence pairs A-B and A-C would be temporarily removed from the analysis ("masked"), leaving only loci that differed between B-C as viable options. Setting maximumPod = 1 prevents this masking step, which can be useful if one is only interested in loci that offer the most discrimination, regardless of what locus pairs offer that discrimination.

(5) Once a locus has been chosen:

a) the specified number of loci has been reached (all unique fingerprints in the case of 'best') and the algorithm terminates; or

b) the specified number of loci has not been reached and there are remaining fingerprints possible, or sequence pairs for which differences exist. The algorithm returns to (1); or

c) there are no remaining fingerprints possible and no sequence pairs for which differences exist. At such time, all sequence pairs are again considered part of the analysis ("unmasked"). If no differences among any sequence pairs exist at this point, the algorithm terminates; if differences remain, the algorithm returns to (1). 

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing@gmail.com)

=cut

package Modules::LociSelector::LociSelector;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Log::Log4perl;

#object creation
sub new {
    my ($class) = shift;
    my $self = {};
    bless( $self, $class );
    $self->_initialize(@_);
    return $self;
}

sub _initialize{
    my($self)=shift;
    my %params=@_;

    #logging
    $self->logger(Log::Log4perl->get_logger());    
    $self->logger->info("Logger initialized in Modules::LociSelector::LociSelector");  

   #on object construction set all parameters
    foreach my $key(keys %params){
        if($self->can($key)){
            $self->$key($params{$key});
        }
        else{
            $self->logger->fatal("$key is not a valid parameter in Modules::LociSelector::LociSelector");
            exit(1);
        }
    }   

    #init data structurs
    $self->_missingChars({});
    $self->_matchedPairs({});
    $self->_data({});
    $self->_selectedLoci({});
    $self->_currentFingerprints([]);
    $self->_dataHeader([]);

    #defaults
    $self->_setMissingChars(['-','?','.']);
    $self->_exhaustedFingerprints(0);
    $self->_numberOfFingerprints(0);
    unless(defined $self->maximizePod){
        $self->maximizePod(0);
    }
}

sub logger{
    my $self=shift;
    $self->{'_logger'}=shift // return $self->{'_logger'};
}

sub _missingChars{
    my $self=shift;
    $self->{'__missingChars'}=shift // return $self->{'__missingChars'};
}

sub _matchedPairs{
    my $self=shift;
    $self->{'__matchedPairs'}=shift // return $self->{'__matchedPairs'};
}


sub inputFile{
    my $self=shift;
    $self->{'_inputFile'}=shift // return $self->{'_inputFile'};
}

sub lociNumber{
    my $self=shift;
    $self->{'_lociNumber'}=shift // return $self->{'_lociNumber'};
}

sub _fileHandle{
    my $self=shift;
    $self->{'__fileHandle'}=shift // return $self->{'__fileHandle'};
}

sub _data{
    my $self=shift;
    $self->{'__data'}=shift // return $self->{'__data'};
}

sub _currentFingerprints{
    my $self=shift;
    $self->{'__currentFingerprints'}=shift // return $self->{'__currentFingerprints'};
}

sub _selectedLoci{
    my $self=shift;
    $self->{'__selectedLoci'}=shift // return $self->{'__selectedLoci'};
}

sub _exhaustedFingerprints{
    my $self=shift;
    $self->{'__exhaustedFingerprints'}=shift // return $self->{'__exhaustedFingerprints'};
}

sub _dataHeader{
    my $self=shift;
    $self->{'__dataHeader'}=shift // return $self->{'__dataHeader'};
}

sub _numberOfFingerprints{
    my $self=shift;
    $self->{'__numberOfFingerprints'}=shift // return $self->{'__numberOfFingerprints'};
}

sub maximizePod{
    my $self=shift;
    $self->{'_maximizePod'}=shift // return $self->{'_maximizePod'};
}

sub run{
    my $self =shift;

    $self->_storeDataInMemory($self->inputFile);

    while(my $locus = $self->_getNextLocus()){
        $self->logger->info("Selected locus $locus");
        $self->_selectedLoci->{$locus}=1;
        $self->_addLocusToCurrentFingerprint($self->_data->{$locus});

        unless($self->maximizePod){
            $self->_updateMatchedPairs($self->_data->{$locus});
        }
    }
    continue{
        if(($self->lociNumber eq 'best') && ($self->_exhaustedFingerprints == 1) ){
            $self->logger->info("Stopping as no new fingerprints are possible. Best number of loci reached.");
            last;
        }
        elsif(scalar(keys %{$self->_selectedLoci}) == $self->lociNumber){
            $self->logger->info("Stopping as user-specified loci number " . $self->lociNumber . " reached");
            last;
        }
    }

    #do the output
    $self->_printResults([sort keys %{$self->_selectedLoci}]);
}

=head2 _updateMatchedPairs

We need to mask the locus pairs that have previously provided discrimination.
This only affects the calculation of the POD for a locus, not the unique fingerprint.
If maximizePod = 1, this step is skipped, and no locus pairs are masked.
Masking the pairs ensures that differences between all columns ("strains") are incorporated.
If one wished to retrieve the loci with the highest POD, regardless of whether locus pairs have
previously been used, the masking can be skipped.

=cut

sub _updateMatchedPairs{
    my $self=shift;
    my $locus=shift;

    my $numberOfLoci=scalar(@{$locus});

    for my $i(0..($numberOfLoci-2)){
        for my $j(($i+1)..($numberOfLoci-1)){           
            unless(defined $self->_missingChars->{$locus->[$i]} || defined $self->_missingChars->{$locus->[$j]}){
                if($locus->[$i] ne $locus->[$j]){
                    $self->_matchedPairs->{$i}->{$j}=1;
                }                              
            }
        }
    }
}

=head2 _addLocusToCurrentFingerprint

Need to add the values from the current locus to the growing fingerprints

=cut

sub _addLocusToCurrentFingerprint{
    my $self=shift;
    my $locus=shift;

    $locus = $self->_substituteMissingChars($locus);
    for my $i(0..(scalar(@{$locus})-1)){
        $self->_currentFingerprints->[$i] = $self->_currentFingerprints->[$i] . $locus->[$i];
    }
}


sub _printResults{
    my $self = shift;
    my $loci = shift;

    $self->logger->info('Printing results');

    foreach my $header(@{$self->_dataHeader}){
        print "\t" . $header;
    }
    print "\n";

    foreach my $locus(@{$loci}){
        print $locus . "\t" . "@{$self->_data->{$locus}}" . "\n";
    }
}

=head2 _setMissingChars

Sets the characters to be ignored when calculating POD and fingerprints.

=cut

sub _setMissingChars{
    my $self=shift;
    my $chars=shift;

    foreach my $char(@{$chars}){
        $self->logger->debug("Setting missing char " . $char);
         $self->_missingChars->{$char}=1;
    }
}


=head2 _getNextLocus

Iterates through the remaining loci until the lociNumber has been hit.

=cut

sub _getNextLocus{
    my $self=shift;

    my $loci;
    unless($self->_exhaustedFingerprints == 1){
        $loci = $self->_getAllBestFingerprintLoci();
    }  

    if(!defined $loci->[0]){
        $self->_exhaustedFingerprints(1);
        $self->logger->info("Fingerprints exhausted");

        if($self->lociNumber eq 'best'){
            return undef;
        }

        $loci = [sort keys %{$self->_data}];
    }
   
    $self->logger->info("Sending " . scalar(@{$loci} . " to _getAllBestPodLoci"));
    my $podLoci = $self->_getAllBestPodLoci($loci);      

    unless(defined $podLoci->[0]){
        $self->logger->info("Resetting matched pairs");
        $self->logger->info("Sending " . scalar(@{$loci} . " to _getAllBestPodLoci"));
        $self->_matchedPairs({});
        $podLoci = $self->_getAllBestPodLoci($loci);
    }
    return $podLoci->[0] // undef;
}


=head2 _getAllBestPodLoci

Given a list of input loci, return all that have the maximum POD,
taking into account those that have already been matched in $self->_matchedPairs

=cut

sub _getAllBestPodLoci{
    my $self=shift;
    my $loci = shift;

    my @bestLoci;
    my $topPod=0;
    foreach my $locus(@{$loci}){

        if(defined $self->_selectedLoci->{$locus}){
            next;
        }    
            
        my $pod = $self->_calculatePod($self->_data->{$locus});
       # $self->logger->debug("Pod locus: $locus pod: $pod");
        if($pod > $topPod){
            @bestLoci=($locus);
            $topPod=$pod;
        }
        elsif($pod == $topPod && $pod !=0){
            push @bestLoci, $locus;
        }
    }
    $self->logger->debug("Top pod: $topPod: numberOfLoci: " . scalar(@bestLoci));
    return \@bestLoci;
}

=head2 _calculatePod

We only need to compare each set once, not bi-directionally.
Thus {i}->{j} is sufficient, we don't need {j}->{i}.
eg. [0][1][2]
The outer loop starts at [0] and goes to [1].
The inner loop starts at [1] and goes to [2].
This gives us:
[0]->[1], [0]->[2]
[1]->[2]
Which is all we need, and saves the needless duplication if both
loops were to run through everything.

=cut

sub _calculatePod{
    my $self=shift;
    my $locus=shift;

    my $pod=0;
    my $numberOfLoci = scalar(@{$locus});
    #$self->logger->debug("Number of pod loci: $numberOfLoci");
    for my $i(0..($numberOfLoci-2)){
        for my $j(($i+1)..($numberOfLoci-1)){
            if(defined $self->_matchedPairs->{$i}->{$j}){
                #$self->logger->debug("Matched pair defined: $i:$j");
                next;
            }
            else{
               # $self->logger->debug("Incrementing pod");
                unless(defined $self->_missingChars->{$locus->[$i]} || defined $self->_missingChars->{$locus->[$j]}){
                    if($locus->[$i] ne $locus->[$j]){
                        $pod++;
                    }
                }                
            }
        }
    }
    #$self->logger->debug("locus: @{$locus} pod: $pod");
    return $pod;  
}

=head2 _getAllBestFingerprintLoci

Looks at all available loci, and returns all that create the most fingerprints from the data.
We want to choose loci that do not contain "missingCharacters" if possible, which is accomplished
by calculating the POD for any set with more than one locus upon return.

=cut

sub _getAllBestFingerprintLoci{
    my $self=shift;  

    $self->logger->info("Getting best fingerprints");

    my @chosenLoci;
    my $initialValue = $self->_numberOfFingerprints;
    my $topValue=$initialValue;

    foreach my $locus(sort keys %{$self->_data}){
        #$self->logger->debug("Next locus: $locus");
        if(defined $self->_selectedLoci->{$locus}){
            next;
        }

        my $fingerprintValue = $self->_calculateFingerprint($self->_data->{$locus});

        if($fingerprintValue > $topValue){
            @chosenLoci=($locus);
            $topValue=$fingerprintValue;
            #$self->logger->debug("New best locus: $locus value: $fingerprintValue");
        }
        elsif($fingerprintValue == $topValue && $topValue != $initialValue){
            push @chosenLoci, $locus;
        }    
    }
    $self->logger->info("Best fingerprint value of $topValue");
    $self->logger->info("initialValue: $initialValue topValue: $topValue");
    $self->_numberOfFingerprints($topValue);
    return \@chosenLoci;
}


=head2 _substituteMissingChars

Determines whether or not a locus contains a missing char.
Replaces with a '.' if present.

=cut

sub _substituteMissingChars{
    my $self=shift;
    my $locus=shift;
  
    for my $i(0..scalar(@{$locus})-1){           
        if(defined $self->_missingChars->{$locus->[$i]}){
            $locus->[$i]='.';
        }
    }
    return $locus;
}

=head2

Returns the number of unique fingerprints, if the current locus is chosen.

=cut

sub _calculateFingerprint{
    my $self = shift;
    my $locus = shift;

    my %tempData;
    foreach my $i(0..(scalar(@{$locus})-1)){
        my $locusValue = $locus->[$i];
        my $value;

        if(defined $self->_missingChars->{$locusValue}){
            $value='.';
            #$self->logger->debug("Missing $value");
        }
        else{
            $value = $locusValue;
        }
        my $key = $self->_currentFingerprints->[$i] . $value;
        $tempData{$key} = 1;
    }
    return $self->_accountForMissingChararcters([keys %tempData]);
}

=head2 _accountForMissingCharacters

We need to ensure that the 'missing characters' are not
seen as unique values, which will give an incorrect, inflated number of
fingerprints.

=cut

sub _accountForMissingChararcters{
    my $self=shift;
    my $prints = shift;

    my $numberOfFingerprints = scalar @{$prints};
    #$self->logger->debug("Accounting for missing chars: number of fingerprints: $numberOfFingerprints");

    for my $i(0..($numberOfFingerprints-2)){
        my $iPrint = $prints->[$i];
        #$self->logger->debug("iprint: $iPrint");
        for my $j(($i+1)..($numberOfFingerprints-1)){
            my $jPrint = $prints->[$j];
            #$self->logger->debug("jPrint: $jPrint");
            if($jPrint =~ m/$iPrint/ || $iPrint =~ m/$jPrint/){
                #$self->logger->debug("Decreasing number of fingerprints from $numberOfFingerprints");
                $numberOfFingerprints--;
            }
            else{
                #$self->logger->debug("$jPrint not equal to $iPrint");
            }          
        }
    }
    return $numberOfFingerprints;
}


=head2 _storeDataInMemory

Stores the input file as a hash, with the locus name as the key,
and the data as an array reference.
Stores the column headers as an arrayref in $self->_dataHeader
for use in output.

=cut

sub _storeDataInMemory{
    my $self=shift;
    my $file =shift;

    my $inFH=IO::File->new('<' . $file) or die "$!";
    $self->logger->info("Collecting data from $file");

    while(my $line = $inFH->getline){ 
        $line =~ s/\R//g;
        my @la = split('\t',$line);

         if($inFH->input_line_number == 1){
            foreach my $head(@la){
                if($head eq ''){
                    next;
                }
                push @{$self->_dataHeader},$head;
            }
            next;
        }

        my $locus = shift(@la);
        $self->_data->{$locus}=\@la;
    }
    $inFH->close();
}

1; #this sentence is not false