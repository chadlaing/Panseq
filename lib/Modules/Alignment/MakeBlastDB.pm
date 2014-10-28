#!/usr/bin/env perl

package Modules::Alignment::MakeBlastDB;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->initialize(@_);
	return $self;
}

#variables
sub programDirectory{
	my $self=shift;
	$self->{'_programDirectory'}=shift // return $self->{'_programDirectory'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub blastDirectory{
	my $self=shift;
	$self->{'_blastDirectory'}=shift // return $self->{'_blastDirectory'};
}

#methods
sub initialize{
	my($self)=shift;
	

	#logging
	$self->logger( Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::MakeBlastDB");

	#init values
	my %params = @_;

	foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			$self->logger->logconfess("$key is not a valid parameter in Modules::Alignment::MakeBlastDB");
		}
	}


	#check required
	unless(
		defined $self->dbtype
		&& defined $self->out
		&& defined $self->title
		&& defined $self->in
		&& defined $self->logfile
		&& defined $self->blastDirectory
	){
		$self->logger->logconfess("Modules::Alignment::MakeBlastDB requires\n
			\t'dbtype'\n
			\t'out'\n
			\t'title'\n
			\t'in'\n
			\t'logfile'\n
			to be defined. One or more parameters missing
		");
	}
}

sub dbtype{
	my $self=shift;
	$self->{'_dbType'}=shift // return $self->{'_dbType'};
}

sub out{
	my $self=shift;
	$self->{'_out'}=shift // return $self->{'_out'};
}

sub title{
	my $self=shift;
	$self->{'_title'}=shift // return $self->{'_title'};
}

sub in{
	my $self=shift;
	$self->{'_in'}=shift // return $self->{'_in'};
}

sub logfile{
	my $self=shift;
	$self->{'_logfile'}=shift // return $self->{'_logfile'};
}


sub run{
	my($self)=shift;
	
	my $systemLine =$self->blastDirectory 
		. 'makeblastdb'
		. ' -dbtype ' . $self->dbtype
		. ' -out '. $self->out
		. ' -max_file_sz 1000000000000' 
		. ' -title ' . $self->title
		. ' -in ' . $self->in
		. ' -logfile ' . $self->logfile;
	
	$self->logger->debug("Running makeblastdb with the following: $systemLine");
	system($systemLine); 
}

1;
