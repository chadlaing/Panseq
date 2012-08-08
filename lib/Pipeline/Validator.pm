#!/usr/local/bin/perl

#Validation Module

#Handles common validation situations. Goal is to make validation more robust.

#Some methods will attempt to compensate for common mistakes or
#alternate but valid inputs by changing it to a standard form.
#eg. adding forward slashes to ends directory names missing the slash

#Constantly added to when new validation situation arise. Also should update
#autocorrection often to make valiation to make processes more robust

package Pipeline::Validator;

use strict;
use diagnostics;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use Carp;
use Scalar::Util qw/ looks_like_number /;

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;
	$self->_initialize(@_);
	return $self;
}

sub _context{
	my $self=shift;
	$self->{'__context'}=shift // return $self->{'__context'};
}

sub _initialize{
	my $self=shift;
	my $caller = caller();
	$self->_context($caller);
}
#setters
sub setContext {
	my $self = shift;
	$self->_context(shift);
	unless ( $self->_context ) {
		croak 'no new context specifed.';
	}
}

#validation methods
sub isDirectoryEmpty {
	my $self=shift;
    my $dir = shift;

    opendir my $FH, $dir or die "Cannot open $dir$!\n";

    while ( defined (my $entry = readdir $FH) ) {
        return unless $entry =~ /^[.][.]?\z/;
    }
    return 1;
}


sub isADirectory {
	my ($self) = shift;

	if (@_) {
		my $dirName = shift;
		if ( $dirName =~ /^\// ) {
			unless ( $dirName =~ /\/$/ ) {

				#add trailing slash if absent
				$dirName .= '/';
			}
			return $dirName;
		}
		else {
			croak $self->_context . ": $dirName directory names must start with '/' in configuration file!\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to isADirectory!\n";
	}
}

sub doesDirectoryExist {
	my ($self) = shift;
	if (@_) {
		my $dirName = shift;
		$dirName = $self->isADirectory($dirName);
		if ( -e $dirName ) {
			return $dirName;
		}
		else {
			croak $self->_context . ": $dirName directory not found\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to doesDirectoryExist!\n";
	}
}

sub isAFilename {
	my ($self) = shift;
	if (@_) {
		my $fileName = shift;
		return $fileName;
	}
	else {
		croak $self->_context . ": nothing sent to isAFilename!\n";
	}
}

sub doesFileExist {
	my ($self) = shift;
	if (@_) {
		my $fileName = shift;
		$fileName = $self->isAFilename($fileName);
		if ( -e $fileName ) {
			return $fileName;
		}
		else {
			croak $self->_context . ": $fileName file not found\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to doesFileExist!\n";
	}
}

sub isAnInt {
	my ($self) = shift;

	if (@_) {
		my $value = shift;

		if ( $value =~ /^-?\d+$/ ) {
			return $value;
		}
		else {
			croak $self->_context . ": $value is not an integer and should be in configuration file! Check for trailing spaces.\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to isAnInt!\n";
	}
}

#Tolerates percentage in decimal or whole number but will convert value to decimal form.
#TODO: In the unlikely case of values of 1, this method will unfortunately take the percentage value
#need to think of clever way to deal with this more elegantly
sub isAValidPercentID {
	my ($self) = shift;

	if (@_) {
		my $number = shift;

		if ( ( $number <= 1 ) && ( $number >= 0 ) ) {
			return $number;
		}

		elsif ( ( $number <= 100 ) ) {
			return $number / 100;
		}
		else {
			croak $self->_context
			  . ":$number is an invalid entry for isAValidPercentID.\n",
			  "Values must be between 0 and 1";

		}
	}
	else {
		croak $self->_context . ": Nothing sent to isAValidPercentID!\n";

	}
}

#accepts yes or no, also 0 or 1 or y and n
#converts input to boolean value (1 or 0)
sub yesOrNo {
	my ($self) = shift;

	if (@_) {
		my $type = shift;
		if ( ( $type eq '1' ) || ( $type eq 'yes' ) || ( $type eq 'y' ) ) {
			return 'yes';
		}
		elsif ( ( $type eq 'no' ) || ( $type eq 'n' ) || ( $type eq '0' ) ) {
			return 'no';
		}
		else {
			croak $self->_context
			  . ": $type is an invalid type for yesOrNo!\n",
			  "Valid types are yes or no, 1 or 0.\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to yesOrNo!\n";

	}
}

#methods taken from panseq Shared and the like

sub accessoryTypeCheck {
	my ($self) = shift;

	if (@_) {
		my $type = shift;

		if ( ( $type eq 'binary' ) || ( $type eq 'percent' ) || ( $type eq 'sequence' ) ) {
			return $type;
		}
		else {
			croak $self->_context . ":$type is not a valid accessoryType!\n", "Valid options are binary, percent and sequence!\n";
		}
	}
	else {
		croak $self->_context . ":nothing sent to accessoryTypeCheck!\n";
	}
}

sub coreComparisonTypeCheck {
	my ($self) = shift;

	if (@_) {
		my $type = shift;

		if ( ( $type eq 'blast' ) || ( $type eq 'nucmer' ) ) {
			return $type;
		}
		else {
			croak $self->_context . ":$type is not a valid coreComparisonType!\n", "Valid options are blast and nucmer!\n";
		}
	}
	else {
		croak $self->_context . ":nothing sent to accessoryTypeCheck!\n";
	}
}

sub blastTypeCheck {
	my ($self) = shift;

	if (@_) {
		my $type = shift;

		if ( $type eq 'blastn' ) {
			return $type;
		}
		elsif ( $type eq 'tblastn' ) {
			return $type;
		}
		else {
			croak $self->_context . ":$type is an invalid entry for blastType in the configuration file!\n
				Currently supported values are blastn and tblastn";
		}
	}
	else {
		croak $self->_context . ":Nothing sent to blastTypeCheck!\n";
	}
}

sub novelRegionFinderModeCheck {
	my ($self) = shift;

	if (@_) {
		my $mode = shift;

		if ( ( $mode eq 'no_duplicates' ) || ( $mode eq 'common_to_all' ) || ( $mode eq 'unique' ) ) {
			return $mode;
		}
		else {
			croak "$mode is not a valid novelRegionFinderMode!\n", "valid modes are no_duplicates, common_to_all and unique.\n";
		}
	}
	else {
		croak $self->_context . ":nothing sent to novelRegionFinderModeCheck\n";
	}
}

#phyML validation methods
#@TODO: *YAGNI Add amino acid models?
sub isEquilibriumModel {
	my ($self) = shift;
	if (@_) {
		my $value = shift;

		if ( $value =~ /[Hh][Kk][Yy]85/ ) {
			return 'HKY85';
		}
		elsif ( $value =~ /[Jj][Cc]69/ ) {
			return 'JC69';
		}
		elsif ( $value =~ /[kK]80/ ) {
			return 'K80';
		}
		elsif ( $value =~ /[Ff]81/ ) {
			return 'F81';
		}
		elsif ( $value =~ /[Ff]84/ ) {
			return 'F84';
		}
		elsif ( $value =~ /[tT][Nn]93/ ) {
			return 'TN93';
		}
		elsif ( $value =~ /[Gg][Tt][Rr]/ ) {
			return 'GTR';
		}
		else {
			croak $self->_context . ": $value not recognized. Equilibrium frequencies must be specified or set to e (emperical) or m (model)!\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to isEqualibrimFreq!\n";
	}
}

sub isPhyMLSampling {
	my ($self) = shift;
	if (@_) {
		my $value = shift;

		if (   ( $value eq 'BEST' )
			|| ( $value eq 'NNI' )
			|| ( $value eq 'SPR' ) )
		{
			return $value;
		}
		else {
			croak $self->_context . ": $value not recognized. Valid Sampling method notations are 'BEST,'NNI' and 'SPR' !\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to isPhyMLSampling!\n";
	}
}

sub equilEstiValue {
	my $self = shift;
	if (@_) {
		my $value = shift;
		if ( $value eq 'e' ) {
			return 'e';
		}
		else {
			return $self->isGreaterThan( $value, 0, 0 );
		}
	}
	else {
		croak $self->_context . ":nothing sent to equilEstiValue\n";
	}
}

sub phyML_v {
	my $self = shift;
	if (@_) {
		my $value = shift;
		if ( $value eq 'e' ) {
			return 'e';
		}
		else {
			return $self->isBetween( $value, 0, 1, 1 );
		}
	}
	else {
		croak $self->_context . ":nothing sent to phyML_v\n";
	}
}

#Also used for phylip
sub isEqualibrimFreq {
	my ($self) = shift;
	if (@_) {
		my $value      = shift;
		my $startValue = shift;

		if (   ( $value eq 'e' )
			|| ( $value eq 'm' )
			|| ( $value =~ /0\.\d+ 0\.\d+ 0\.\d+ 0\.\d+/ ) )
		{
			return $value;
		}
		else {
			croak $self->_context . ": $value not recognized. Equilibrium frequencies must be specified or set to e (emperical) or m (model)!\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to isEqualibrimFreq!\n";
	}
}

#mrBayes unique validations

sub isRateVarType {
	my ($self) = shift;
	if (@_) {
		my $value = shift;

		if (   ( $value eq 'equal' )
			|| ( $value eq 'gamma' )
			|| ( $value eq 'propinv' )
			|| ( $value eq 'invgamma' )
			|| ( $value eq 'adgamma' ) )
		{
			return $value;
		}
		else {
			croak $self->_context
			  . ": $value not recognized. Rate variation must be specified or set to equal, gamma, propinv, invgamma or adgamma!\n";
		}
	}
	else {
		croak $self->_context . ": nothing sent to isRateVarType!\n";
	}
}

#Advanced methods (these require more than one input)
#requires a secondary number to validate greater than.
#usage isGreaterThan($value, $startValue, $isInclusive);
sub isGreaterThan {
	my ($self) = shift;

	if (@_) {
		my $value       = shift;
		my $startValue  = shift;
		my $isInclusive = shift;

		if ($isInclusive) {
			if ( looks_like_number($value)
				&& ( $value >= $startValue ) )
			{
				return $value;
			}
			else {
				croak $self->_context . ": $value should be a value greater than $startValue in configuration file!\n";
			}
		}
		else {
			if ( looks_like_number($value)
				&& ( $value > $startValue ) )
			{
				return $value;
			}
			else {
				croak $self->_context . ": $value should be a value greater than $startValue in configuration file!\n";
			}
		}
	}
	else {
		croak $self->_context . ": nothing sent to isGreaterThan!\n";
	}
}

#requires a secondary number to validate greater than.
sub isAnIntGreaterThan {
	my ($self) = shift;
	if (@_) {
		my $value      = shift;
		my $startValue = shift;
		$value = $self->isAnInt($value);
		$value = $self->isGreaterThan( $value, $startValue );
		return $value;
	}
	else {
		croak $self->_context . ": nothing sent to isAnIntGreaterThan!\n";
	}
}

#requires a secondary lower bound value and a tertiary upper bound value
#usage isBetween($value, $startValue, $endValue,$isInclusive);
sub isBetween {
	my ($self) = shift;
	if (@_) {
		my $value       = shift;
		my $startValue  = shift;
		my $endValue    = shift;
		my $isInclusive = shift;
		if ($isInclusive) {
			if (   looks_like_number($value)
				&& ( $value >= $startValue )
				&& ( $value <= $endValue ) )
			{
				return $value;
			}
			else {
				croak $self->_context . ": $value should be a value between $endValue and $startValue in configuration file!\n";
			}
		}
		else {
			if (   looks_like_number($value)
				&& ( $value > $startValue )
				&& ( $value < $endValue ) )
			{
				return $value;
			}
			else {
				croak $self->_context . ": $value should be a value between $endValue and $startValue in configuration file!\n";
			}
		}
	}
	else {
		croak $self->_context . ": nothing sent to isBetween!\n";
	}
}

#requires a secondary lower bound value and a tertiary upper bound value
sub isAnIntBetween {
	my ($self) = shift;
	if (@_) {
		my $value      = shift;
		my $startValue = shift;
		my $endValue   = shift;
		$value = $self->isAnInt($value);
		$value = $self->isBetween( $value, $startValue, $endValue );
		return $value;
	}
	else {
		croak $self->_context . ": nothing sent to isAnIntBetween!\n";
	}
}

sub isAValidWebpageType{
	my $self=shift;
	my $type = shift // 'empty';
	
	if(
		($type eq 'home') ||
		($type eq 'novel') ||
		($type eq 'core') ||
		($type eq 'loci') ||
		($type eq 'tutorial') ||
		($type eq 'contact') ||
		($type eq 'submit')
		
	){
		return $type
	}
	else{
		croak $self->_context . ": incorrect web page type!\n";
	}
}

1;
