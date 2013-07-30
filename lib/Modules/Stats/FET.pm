#!/usr/bin/env perl
package Modules::Stats::FET;

=pod

=head1 NAME

Stats::FET - Given a tab-delimited table of data, and a list of two groups, returns FET p-value of target character between the two groups for every row of the table.

=head1 SYNOPSIS

	my $fet = Stats::FET->new(
		'testCharacters'=>['A','C','T','G'],
		'inputFile'=>'/your/file/name.txt',
		'excludedCharacters'=>['-'],
		'group1'=>['strain1','strain3','strain5','strain6','strain88'],
		'group2'=>['strain131','strain120','strain111','strain44','strain55','strain66','strain77']
	);

	OR

	my $fet = Stats::FET->new(
		'testCharacters'=>['A','C','T','G'],
		'inputFile'=>'/your/file/name.txt',
		'excludedCharacters'=>['-'],
		'group1'=>'/file/listing/group1.txt',
		'group2'=>'/file/listing/group2.txt'
	);

	$fet->run();


=head1 DESCRIPTION

Stats::FET is for analyzing tab-delimited table data with both column and row labels.
An object is created by specifying the file name of the table one will analyze, along with single character data to test for significance with,
characters that should be ignored rather than counted as "not the test char", and the names of the members of two groups.
The column headers must match exactly the group names specified upon object creation.
The module prints the results sorted by ascending P-value as follows:

<row header>	<character tested for significance>	<P-value>

The default output is STDOUT, but can be piped or specified directly using $obj->outputFH(<filehandle>);

=head2 Methods

=head3 inputFile

Specifies the tab-delimited table data file to be analyzed.
Must be specified during object creation.
RO.

=head3 _setInputFile

Private method to set inputFile at object creation.

=head3 testCharacters

An ARRAYREF to the single characters to be tested for significance between the two groups.
Must be specified during object creation.
RO.

=head3 _setTestCharacters

Private method to set testCharacters at object creation.

=head3 group1

A HASHREF that stores the members of group1.
$self->group1->{strain}=1/
Must be specified during object creation.
RO.

=head3 _group1Columns

Private method that stores each column number for strains in group 1 in an ARRAYREF.

=head3 _group2Columns

Private method that stores each column number for strains in group 2 in an ARRAYREF.

=head3 group2

A HASHREF that stores the members of group2.
$self->group2->{strain}=1/
Must be specified during object creation.
RO.

=head3 _resultHash

Private HASHREF that stores the results of the entire run.
$self->_resultHash->{row_header}->{test_character}=p_value.

=head3 _memoize

A HASHREF that stores the four number FET query sequence as the hash key and P Value for each FET as the value.
	$self->_memoize->{'2_4_8_11'}=p_value
This allows the program to check if a test has previously been computed and if so fetch the P Value from the hash.
Both combinations that give the same P Value are automatically stored, eg.
FET(2,4,8,11) == FET(8,11,2,4)

=head3 logger

Stores the Log::Log4perl object for module logging.

=head3 _R

Private method that stores the Statistics::R object and is used for all R-related calls.

=head3 _excludedCharactersHash

Private HASHREF to all characters that should not be included in the FET, rather than being counted as !testCharacter.

=head3 _initialize

Private method called by new to initialize object.
Calls the FileInteraction::FlexiblePrinter _initialize.
Sets the Log::Log4perl object.
Sets the Statistics::R object.
Sets all required class variables:
	inputFile
	testCharacters
	excludedCharacters
	group1
	group2

Initializes the data structures:
	$self->_memoize({});
	$self->_resultHash({});

=head3 _setExcludedCharacters

Private method to set excludedCharacters at object creation.

=head3 _setGroup

	$self->_setGroup('grpup_name',(ARRAYREF or <filename>);

Private method that converts an ARRAYREF or items listed in a file into the _group1 or _group2 HASHREF, with the key as the group member name and the value of 1.

=head3 _getPValue

	$self->_getPValue(group1pos,group1neg,group2pos,group2neg);
Private method that returns the P value for the 4 numbers passed as an ARRAY to the function.
Uses the two-tailed test via $self->_R and fisher.test().

=head3 _setGroupColumns

Private method that reads the first row of the table, and starting at the second column (to allow the row headers to occupy the first column)
creates the ARRAYREF _group1Columns or _group2Columns that have the column numbers for each group member stored.

=head3 _processLine

	$self->_processLine($line,$testChar);

Private method that prepares the four values to be sent to the _getPValue function.
Retrieves the counts for each group via the _countGroup function.
Queries _memoize to see if FET has previously been computed, if it has, gets P Value from hash;
if not, calls the _getPValue function and adds the result to the _memoize HASHREF.
Stores results in the _resultHash.

=head3 _countGroup

	$self->_countGroup('groupname', $testChar, \@la);

Where @la is the array of the line, with each array cell containing a single character, computed by the _processLine function where this function is called.
Returns the hashRef where $hashRef->{'pos'}=(number of times the testCharacter was counted)
and $hashRef->{'neg'}=(number of non-testCharacters, excluding any _excludedCharacters).

=head3 _printResultsToFile

Sorts the _resultHash by ascending P Value and outputs the results in the form:

row_header	test_character	count_data	P_Value

uses the inherited ->print and ->outputFH methods from FlexiblePrinter, which defaults to STDOUT.

=head3 run

The public method to get the FET stats from the file.
Opens the file, iterates over each line, and send the line to the _processLine function, except for the first line which defines the column arrays via _setGroupColumns.
Finally, calls the _printResultsToFile function.

=head3 DESTROY

Closes the _R instance started in _initialize.


=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing gmail com)

=cut
use FindBin;
use lib "$FindBin::Bin../";
use warnings;
use strict;
use Carp;
use IO::File;
use Statistics::R;
use Log::Log4perl;
use Role::Tiny::With;

with 'Roles::FlexiblePrinter';

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

sub inputFile{
	my $self=shift;
	return $self->{'_inputFile'};
}

sub _setInputFile{
	my $self=shift;
	$self->{'_inputFile'}=shift;
}


sub testCharacters{
	my $self=shift;
	return $self->{'_testCharacters'};
}

sub _setTestCharacters{
	my $self=shift;
	$self->{'_testCharacters'}=shift;
}

sub group1{
	my $self=shift;
	return $self->{'_group1'};
}

sub _group1Columns{
	my $self=shift;
	$self->{'__group1Column'}=shift // return $self->{'__group1Column'};
}

sub _group2Columns{
	my $self=shift;
	$self->{'__group2Column'}=shift // return $self->{'__group2Column'};
}

sub group2{
	my $self=shift;
	return $self->{'_group2'};
}

sub _resultHash{
	my $self=shift;
	$self->{'__resultHash'}=shift // return $self->{'__resultHash'};
}

sub _memoize{
	my $self=shift;
	$self->{'__memoize'}=shift // return $self->{'__memoize'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub _R{
	my $self=shift;
	$self->{'__R'}=shift // return $self->{'__R'};
}

sub _excludedCharactersHash{
	my $self=shift;
	$self->{'__excludedCharactersHash'}=shift // return $self->{'__excludedCharactersHash'};
}


#methods
sub _initialize{
	my $self=shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Stats::FET");

	my %init = @_;
	
	$self->_setInputFile($init{'inputFile'}) // confess("inputFile required in FET object creation");
	$self->_setTestCharacters($init{'testCharacters'}) // confess("testCharacters required in FET object creation");
	$self->_setExcludedCharacters($init{'excludedCharacters'})	// confess("excludedCharacters required in FET object creation");
	$self->_setGroup('group1',$init{'group1'}) // confess("group1 required in FET object creation");
	$self->_setGroup('group2',$init{'group2'}) // confess("group2 required in FET object creation");

	#set R
	$self->_R(Statistics::R->new());
	$self->_R->startR;

	#initialize data structures
	$self->_memoize({});
	$self->_resultHash({});

	return 1;
}

sub _setExcludedCharacters{
	my $self=shift;
	my $excludedCharacters=shift;

	my %excludedCharactersHash = map {$_ => 1} @{$excludedCharacters};
	$self->_excludedCharactersHash(\%excludedCharactersHash);

	return 1;
}

sub _setGroup{
	my $self=shift;

	#allow for the input to be either a filename of single words, each on a newline, or an array ref to group names
	#generate a hash of names for quick checking

	my $group=shift;
	my $input = shift;
	my $type = ref($input);
	my $outputHashRef;

	if($type =~ /ARRAY/ ){
		#ARRAYREF
		#groups already defined in ARRAYREF, do nothing to input
		my %strainsHash = map{$_ => 1} @{$input};
		$outputHashRef = \%strainsHash;
	}
	elsif(-e $input){
		#FILE name input
		#create an array and set input as arrayref
		my %strainsHash;

		my $inFH = IO::File->new('<' . $input) or die "cannot open $input $!\n";
		while (my $line = $inFH->getline) {
			$line =~ s/\s//g;
			$strainsHash{$line}=1;
		}
		$inFH->close();
		$outputHashRef = \%strainsHash;
	}
	else{
		$self->logger->fatal("Group file does not exist or incorrect group file");
		exit(1);
	}

	if($group eq 'group1'){
		$self->{'_group1'}=$outputHashRef;
		$self->logger->debug("group1 set to $outputHashRef");
	}
	elsif($group eq 'group2'){
		$self->{'_group2'}=$outputHashRef;
		$self->logger->debug("group2 set to $outputHashRef");
	}
	else{
		$self->logger->error("Incorrect group type (other than group1 or group2) specified in FET");
		exit(1);
	}

	return 1;
}

sub _getPValue{
	my $self=shift;
	my @values=@_;

	my $matrixString =
		    'comp<-matrix(c('
		  . $values[0] . ','
		  . $values[1] . ','
		  . $values[2] . ','
		  . $values[3] . ' ), nr = 2)';

	my $rQuery = 'fisher.test(' . $matrixString . ")\n";
	
	#access R	
	$self->_R->send($rQuery);
	my $results = $self->_R->read();
	
	my $pvalue;
	if($results =~ /p-value\D+(.+)/){		
		$pvalue=$1;
		$self->logger->debug("p-value: $pvalue");
	}
	else{
		$self->logger->fatal("no p-value $results");
		exit(1);
	}
	return $pvalue;
}


sub _setGroupColumns{
	my $self = shift;
	my $line = shift;

	$line =~ s/\R//g;
	my @la = split(/\t/,$line);
	$self->logger->debug("line: $line");

	my @group1Columns;
	my @group2Columns;
	for my $column(1..(scalar(@la)-1)){
		my $name = $la[$column];

		#if there are non-printing chars, check and skip
		$name =~ s/\s//g;
		if((!defined $name) || ($name eq '')){
			next;
		}

		$self->logger->debug("name: $name");

		#add the column numbers for each group to a private array
		if(defined $self->group1->{$name}){
			push @group1Columns, $column;
		}
		elsif(defined $self->group2->{$name}){
			push @group2Columns, $column;
		}
	}
	$self->_group1Columns(\@group1Columns);
	$self->_group2Columns(\@group2Columns);
	return 1;
}

sub _countGroup{
	my $self=shift;
	my $group=shift;
	my $testChar=shift;
	my $laRef = shift;

	#get either group1 or group2
	my @columns;
	if($group eq 'group1'){
		@columns = @{$self->_group1Columns};
	}
	elsif($group eq 'group2'){
		@columns = @{$self->_group2Columns};
	}
	else{
		$self->logger->fatal("Incorrect group sent to _countGroup");
		exit(1);
	}


	#do the counts
	my %countHash;
	$countHash{'pos'}=0;
	$countHash{'neg'}=0;

	foreach my $column(@columns){
		my $locus=$laRef->[$column];

		if(defined $self->_excludedCharactersHash->{$locus}){
			next;
			$self->logger->debug("Skipping locus $locus");
			#this prevents any addition to pos/neg count
		}

		if($locus eq $testChar){
			$countHash{'pos'}++;
		}
		else{
			$countHash{'neg'}++;
		}
	}

	#return the counts in a hash
	return \%countHash; 
}

sub _processLine{
	my $self=shift;
	my $line = shift;
	my $testChar = shift;

	$line =~ s/\R//g;
	my @la = split(/\t/,$line);
	#$self->logger->debug('number of array elements:' . scalar(@la));

	my $locusName = $la[0];
	
	#get the count of all in each group that have the test char
	#check excludedCharacters hash to ensure absent data is not treated the same as !char
	my $group1Counts = $self->_countGroup('group1',$testChar,\@la);
	my $group2Counts = $self->_countGroup('group2',$testChar,\@la);

	#run the FET using the current values if not in _memoize
	my $memoizeString = join('_', $group1Counts->{'pos'}, $group1Counts->{'neg'}, $group2Counts->{'pos'}, $group2Counts->{'neg'});

	#$self->logger->debug('memoizeString: ' . $memoizeString);

	my $pValue;	
	if(defined $self->_memoize->{$memoizeString}){
		$pValue = $self->_memoize->{$memoizeString};
	}
	else{
		#get pValue
		#add to memoize
		#add reverse to memoize as it is the same eg fet(1,4,8,11) == fet(8,11,1,4)
		$pValue = $self->_getPValue($group1Counts->{'pos'}, $group1Counts->{'neg'}, $group2Counts->{'pos'}, $group2Counts->{'neg'});
		$self->_memoize->{$memoizeString}=$pValue;

		my $memoizeString2 = join('_', $group2Counts->{'pos'}, $group2Counts->{'neg'}, $group1Counts->{'pos'}, $group1Counts->{'neg'});
		$self->_memoize->{$memoizeString2} = $pValue;
	}

	#add to result hash
	$self->_resultHash->{$locusName}->{$testChar}=$pValue;

	return 1;
}

sub _printResultsToFile{
	my $self=shift;

	my %outputHash;
	foreach my $locus(keys %{$self->_resultHash}){
		foreach my $testChar(keys %{$self->_resultHash->{$locus}}){
			$outputHash{$locus . "\t" . $testChar} = $self->_resultHash->{$locus}->{$testChar};
		}
	 }

	 #sort new hash by value but output the key
	 foreach my $output(sort{$outputHash{$a} <=> $outputHash{$b}} keys %outputHash){
	 	$self->printOut($output . "\t" . $outputHash{$output} . "\n");
	 }

	 return 1;
}

sub run{
	my $self=shift;

	my $inFH = IO::File->new('<' . $self->inputFile) or die 'Cannot open ' . $self->inputFile . "$!\n";
	while(my $line = $inFH->getline){
		if($.==1){
			$self->_setGroupColumns($line);
			next;
		}
		else{
			#process each line for the FET p-value of the specified groups
			foreach my $testChar(@{$self->testCharacters}){
				#_processLine adds the pvalue to $self->_resultHash->{locus_name}->{test_char}=p_value
				$self->_processLine($line,$testChar);				
			}
		}
	}
	$inFH->close();

	$self->_printResultsToFile();

	return 1;
}


#define a DESTROY method to close R instance on object close
sub DESTORY{
	my $self=shift;
	$self->_R->stopR;
}

1;