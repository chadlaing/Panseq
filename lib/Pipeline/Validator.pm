#!/usr/local/bin/perl

#Validation Module 19/1/12 - JC

#Handles common validation situations. Goal is to make validation more robust.

#Some methods will attempt to compensate for common mistakes or
#alternate but valid inputs by changing it to a standard form.
#eg. adding forward slashes to ends directory names missing the slash

#Constantly added to when new validation situation arise. Also should update
#autocorrection often to make valiation to make processes more robust

package Validator;

use strict;
use diagnostics;
use warnings;
use Carp;
use Object::Tiny::RW ('_context');
use Scalar::Util qw/ looks_like_number /;

#constuctor: invoke with string with context in denoting process or module
#eg. Validation->new('Mummer Core Accessory Analysis');
sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;
	$self->_context(shift);
	unless ( $self->_context ) {
		$self->_context('Generic Validation');
	}
	return $self;
}

#setters
sub setContext {
	my $self = shift;
	$self->_context(shift);
	unless ( $self->_context ) {
		confess 'no new context specifed.';
	}
}

#validation methods
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
			confess $self->_context
			  . ": $dirName directory names must start with '/' in configuration file!\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to isADirectory!\n";
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
			confess $self->_context . ": $dirName directory not found\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to doesDirectoryExist!\n";
	}
}

sub isAFilename {
	my ($self) = shift;
	if (@_) {
		my $fileName = shift;
		if ( $fileName =~ /^\// ) {
			return $fileName;
		}
		else {
			confess $self->_context
			  . ": $fileName filenames must start with '/' in configuration file!\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to isAFilename!\n";
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
			confess $self->_context . ": $fileName file not found\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to doesFileExist!\n";
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
			confess $self->_context
			  . ": $value is not an integer and should be in configuration file! Check for trailing spaces.\n";

		}
	}
	else {
		confess $self->_context . ": nothing sent to isAnInt!\n";
	}
}

#tolerates percentage in decimal or whole number but will convert value to decimal form.
sub isAValidPercentID {
	my ($self) = shift;

	if (@_) {
		my $number = shift;

		if ( ( $number <= 1 ) && ( $number >= 0 ) ) {
			return $number;
		}

		elsif ( ( $number <= 100 ) && ( $number >= 0 ) ) {
			return $number / 100;
		}
		else {
			confess $self->_context
			  . ":$number is an invalid entry for isAValidPercentID.\n",
			  "Values must be between 0 and 1";

		}
	}
	else {
		confess $self->_context . ": Nothing sent to isAValidPercentID!\n";

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
			confess $self->_context
			  . ": $type is an invalid type for yesOrNo!\n",
			  "Valid types are yes or no, 1 or 0.\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to yesOrNo!\n";

	}
}

#checks for acceptible tree generation program/method
sub isTreeGenProgram {
	my ($self) = shift;
	if (@_) {
		my $value = shift;
		if ( $value =~ /[pP][hH][yY][mM][lL]/ ) {
			return 'phyML';
		}
		elsif ( $value =~ /[mM][rR][Bb][aA][yY][eE][sS]/ ) {
			return 'mrBayes';
		}
		elsif ( $value =~ /[pP][hH][yY][mM][lL][nN][jJ]/ ) {
			return 'phylipNJ';
		}
		elsif (( $value =~ /[pP][hH][yY][mM][lL][mM][pP]/ )
			|| ( $value =~ /[pP][hH][yY][mM][pP][aA][rR][sS]/ ) )
		{
			return 'phylipMP';
		}
		elsif ( $value =~ /[pP][hH][yY][mM][lL][bB][pP][nN][jJ]/ ) {
			return 'phylipBPNJ';
		}
		elsif (( $value =~ /[pP][hH][yY][mM][lL][bB][pP][mM][pP]/ )
			|| ( $value =~ /[pP][hH][yY][mM][lL][bB][pP][aA][rR][sS]/ ) )
		{
			return 'phylipBPMP';
		}
		else {
			confess $self->_context
			  . ": $value not recognized as a known Phylogenetic Tree Program/Method!\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to isTreeGenProgram!\n";
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
			confess $self->_context
			  . ": $value not recognized. Equilibrium frequencies must be specified or set to e (emperical) or m (model)!\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to isEqualibrimFreq!\n";
	}
}

sub isPhyMLSampling {
	my ($self) = shift;
	if (@_) {
		my $value = shift;

		if (   ( $self->_s eq 'BEST' )
			|| ( $self->_s eq 'NNI' )
			|| ( $self->_s eq 'SPR' ) )
		{
			return $value;
		}
		else {
			confess $self->_context
			  . ": $value not recognized. Valid Sampling method notations are 'BEST,'NNI' and 'SPR' !\n"
			  ;
		}
	}
	else {
		confess $self->_context . ": nothing sent to isPhyMLSampling!\n";
	}
}

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
			confess $self->_context
			  . ": $value not recognized. Equilibrium frequencies must be specified or set to e (emperical) or m (model)!\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to isEqualibrimFreq!\n";
	}
}

#Advanced methods (these require more than one input)
#requires a secondary number to validate greater than.
sub isGreaterThan {
	my ($self) = shift;

	if (@_) {
		my $value      = shift;
		my $startValue = shift;

		if ( looks_like_number($value)
			&& ( $value > $startValue ) )
		{
			return $value;
		}
		else {
			confess $self->_context
			  . ": $value should be a value greater than $startValue in configuration file!\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to isGreaterThan!\n";
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
	}
	else {
		confess $self->_context . ": nothing sent to isAnIntGreaterThan!\n";
	}
}

#requires a secondary lower bound value and a tertiary upper bound value
sub isBetween {
	my ($self) = shift;
	if (@_) {
		my $value      = shift;
		my $startValue = shift;
		my $endValue   = shift;
		if (   looks_like_number($value)
			&& ( $value > $startValue )
			&& ( $value < $endValue ) )
		{
			return $value;
		}
		else {
			confess $self->_context
			  . ": $value should be a value between $endValue and $startValue in configuration file!\n";
		}
	}
	else {
		confess $self->_context . ": nothing sent to isBetween!\n";
	}
}
1;
