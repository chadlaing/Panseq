#!/usr/bin/perl

package Logger;

use strict;
use warnings;
use diagnostics;
use Log::Log4perl qw(:easy);

#object creation
sub new{
	my($class)=shift;
    my $self = {};
    bless ($self, $class);
    $self->_loggerInitialize(@_);
    return $self;
}

#class variables
sub logger{
	my $self=shift;
	$self->{'_Logger_logger'}=shift // return $self->{'_Logger_logger'};
}

sub _loggerInitialize{
	my $self=shift;
	
	Log::Log4perl->easy_init(
	    {
		    level => $DEBUG,  # one of DEBUG, INFO, WARN, ERROR, FATAL
		    #file  => 'STDERR',
		    file=>'/home/phac/panseq/RUNLOG.txt'
	    }
	);
	my $logger= get_logger();
	$self->logger($logger);	
}

1;