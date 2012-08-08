#!/usr/bin/perl

#10-31-11 - Script intended for use with Muscle.
#Functionality added:
#-in
#-out
#-stable
#> (out pipes)
#< (in pipes)
#-maxiters
#-diags
#-refine

package Muscle::MuscleCmd;

use warnings;
use strict;
use diagnostics;
use Carp;
use FindBin;
use lib "$FindBin::Bin";

#constructor
sub new {
	my ($class) = shift;	
	my $self = {};
	bless $self, $class;
	$self->_initialize(@_);
	return $self;
}

sub _initialize{
	my $self=shift;
	my $executable=shift;
	$self->_muscleExecutable($executable) // confess("The location of the muscle executable must be specified");
	$self->_stable(0);
	$self->_diagsNum(0);
	$self->_refine(0);
}

sub _muscleExecutable{
	my $self=shift;
	$self->{'__muscleExecutable'}=shift // return $self->{'__muscleExecutable'};
}

sub _refine{
	my $self=shift;
	$self->{'__refine'}=shift // return $self->{'__refine'};
}

sub _diagsNum{
	my $self=shift;
	$self->{'__diagsNum'}=shift // return $self->{'__diagsNum'};
}

sub _maxiters{
	my $self=shift;
	$self->{'__maxiters'}=shift // return $self->{'__maxiters'};
}

sub _outPipe{
	my $self=shift;
	$self->{'__outPipe'}=shift // return $self->{'__outPipe'};
}

sub _inPipe{
	my $self=shift;
	$self->{'__inPipe'}=shift // return $self->{'__inPipe'};
}

sub _stable{
	my $self=shift;
	$self->{'__stable'}=shift // return $self->{'__stable'};
}

sub _out{
	my $self=shift;
	$self->{'__out'}=shift // return $self->{'__out'};
}

sub _in{
	my $self=shift;
	$self->{'__in'}=shift // return $self->{'__in'};
}

#set functions
#include file path
sub setIn {
	my $self = shift;
	if ( $self->_inPipe ) { confess 'Piped type in command already used.'; }
	else {
		$self->_in(shift);
		unless ( $self->_in ) {
			confess 'Require input, Usage: setIn($self->_in_filename)';
		}
	}
}

#include file path
sub setOut {
	my $self = shift;
	if ( $self->_outPipe ) { confess 'Piped type out command already used.'; }
	else {
		$self->_out(shift);
		unless ( $self->_out ) {
			confess 'Require input, Usage: setOut($self->_out_filename)';
		}
	}
}

sub setInPipe {
	my $self = shift;
	if ( $self->_in ) { confess 'In command already used.'; }
	else {
		$self->_inPipe(shift);
		unless ( $self->_inPipe ) {
			confess 'Require input, Usage: setInPipe($self->_in_filename)';
		}
	}
}

sub setOutPipe {
	my $self = shift;
	if ( $self->_out ) { confess 'Out command already used.'; }
	else {
		$self->_outPipe(shift);
		unless ( $self->_outPipe ) {
			confess 'Require input, Usage: setOutPipe($self->_out_filename)';
		}
	}
}

sub setDiagsOn {
	my $self = shift;
	$self->_diagsNum(1);
}

sub setMaxiters {
	my $self = shift;
	$self->_maxiters(shift);
	unless ( $self->_maxiters ) {
		confess 'Require input, Usage: setmaxiters($self->_maxiters)';
	}
}

sub setRefineOn {
	my $self = shift;
	$self->_refine(1);
}

#turns on ordering of sequence by order input and not similarity (default off)
sub setStableOn {
	my $self = shift;
	$self->_stable(1);
}

#get functions
sub getIn {
	my $self = shift;
	return $self->_in;
}

sub getOut {
	my $self = shift;
	return $self->_out;
}

sub getOutPipe {
	my $self = shift;
	return $self->_outPipe;
}

sub getInPipe {
	my $self = shift;
	return $self->_inPipe;
}

sub isStable {
	my $self = shift;
	return $self->_stable;
}

sub isDiagsOn {
	my $self = shift;
	return $self->_diagsNum;
}

sub getMaxiters {
	my $self = shift;
	return $self->_maxiters;
}

sub isRefineOn {
	my $self = shift;
	return $self->_refine;
}


#primary method
sub run {
	my $self       = shift;
	my $systemLine = $self->_muscleExecutable . ' -quiet';
	unless ( $self->_in || $self->_inPipe ) {
		confess 'Require input using setIn() or setInPipe()';
	}
	if ( $self->_in ) {
		$systemLine .= ' -in \'' . $self->_in . '\'';
	}
	if ( $self->_out ) {
		$systemLine .= ' -out \'' . $self->_out . '\'';
	}
	if ( $self->_inPipe ) {
		$systemLine .= ' < \'' . $self->_inPipe . '\'';
	}
	if ( $self->_outPipe ) {
		$systemLine .= ' > \'' . $self->_outPipe . '\'';
	}
	if ( $self->_stable == 1 ) {
		$systemLine .= ' -stable';
	}
	if ( $self->_maxiters ) {
		$systemLine .= ' -maxiters ' . $self->_maxiters;
	}
	if ( $self->_diagsNum == 1 ) {
		$systemLine .= ' -diags';
	}
	if ( $self->_refine ) {
		$systemLine .= ' -refine';
	}
	system($systemLine);
}
1;
