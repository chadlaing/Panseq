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

package MuscleCmd;
use warnings;
use strict;
use diagnostics;
use Object::Tiny::RW qw{
  _in
  _out
  _stable
  _inPipe
  _outPipe
  _maxiters
  _diagsNum
  _refine
};

use Carp;

#constructor
sub new {
	my ($class) = shift;
	my $self = {};
	bless $self, $class;

	$self->_stable(0);
	$self->_diagsNum(0);
	$self->_refine(0);
	return $self;
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
	my $systemLine = 'muscle -quiet';
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
