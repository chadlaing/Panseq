#!/usr/bin/env perl

package Modules::Alignment::NucmerRun;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Log::Log4perl;
use Carp;

#object creation
sub new{
	my $class=shift;
    my $self = {};
    bless ($self, $class);
    $self->_initialize(@_);   
    return $self;
}

sub _initialize{
	my $self=shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::NucmerRun\n");

	#init values
	my %params = @_;

	foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			$self->logger->logconfess("$key is not a valid parameter in Modules::Alignment::NucmerRun");
		}	
	}
}

#class variables
sub logFile{
	my $self=shift;
	$self->{'_mummerLog'} = shift // return $self->{'_mummerLog'};
}

sub percentIdentityCutoff{
	my $self=shift;
	$self->{'_percentIdentityCutoff'} = shift // return $self->{'_percentIdentityCutoff'};
}

sub mummerDirectory{
	my $self=shift;
	$self->{'_mummerDirectory'} = shift // return $self->{'_mummerDirectory'};
}

sub queryFile{
	my $self=shift;
	$self->{'_queryFile'}=shift // return $self->{'_queryFile'};
}

sub referenceFile{
	my $self=shift;
	$self->{'_referenceFile'}=shift // return $self->{'_referenceFile'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub numberOfCores{
	my $self=shift;
	$self->{'_numberOfCores'}=shift // return $self->{'_numberOfCores'};
}

sub referenceFileArray{
	my $self=shift;
	$self->{'_referenceFile'}=shift // return $self->{'_referenceFile'};
}

sub refFileSizeLimit{
	my $self=shift;
	$self->{'_refFileSizeLimit'}=shift // return $self->{'_refFileSizeLimit'};	
}

sub systemLineBase{
	my $self=shift;
	$self->{'_systemLineBase'}=shift // return $self->{'_systemLineBase'};
}


sub l{
	my $self=shift;
	$self->{'_l'}=shift // return $self->{'_l'};
}

sub g{
	my $self=shift;
	$self->{'_g'}=shift // return $self->{'_g'};
}

sub d{
	my $self=shift;
	$self->{'_d'}=shift // return $self->{'_d'};
}

sub c{
	my $self=shift;
	$self->{'_c'}=shift // return $self->{'_c'};
}

sub b{
	my $self=shift;
	$self->{'_b'}=shift // return $self->{'_b'};
}

sub p{
	my $self=shift;
	$self->{'_p'}=shift // return $self->{'_p'};
}

sub coordsFile{
	my $self=shift;
	$self->{'_coordsFile'}=shift // return $self->{'_coordsFile'};
}

sub _tempFiles{
	my $self=shift;
	$self->{'__tempFiles'}=shift // return $self->{'__tempFiles'};
}



sub run{
	my $self=shift;
	
	#check for required files
	unless(defined $self->referenceFile){
		$self->logger->logconfess("Reference file required in run\n");
	}

	unless(defined $self->queryFile){
		$self->logger->logconfess("Query file required in run\n");
	}

	unless(defined $self->p){
		$self->logger->logconfess("File prefix required in run\n");
	}

	unless(defined $self->mummerDirectory){
		$self->logger->logconfess("Mummer directory required in run");
	}

	#run mummer
	my $deltaFile = $self->p . '.delta';
	my $mummerLine = $self->_createMummerLine();
	$self->logger->debug("Running nucmer comparison with the command: $mummerLine");
	system($mummerLine);	

	#add a delta-filter to limit the nucmer matches to the percentIdentity cutoff
	#this differs slightly from how BLAST calculates things
	#letting all matches go through, to remove discrepancy
	#my $filteredDeltaFile = $self->_deltaFilter($deltaFile);
	
	#run show-coords on delta file
	$self->_showCoords($deltaFile);
}

=head2

Filter the delta-file to contain only those matches that are at or
above the percentIdentityCutoff.

=cut

sub _deltaFilter{
	my $self=shift;
	my $deltaFile = shift;

	#delta-filter [options] <delta file> > <filtered delta file>
	my $filteredDeltaFile = $self->p . '_filtered.delta';
	my $filterLine = $self->mummerDirectory . 'delta-filter -i ' . $self->percentIdentityCutoff
		. ' ' . $deltaFile . ' > ' . $filteredDeltaFile;

	$self->logger->debug("Filtering delta-file for minimum percent identity of " . $self->percentIdentityCutoff
	 . " with command:\n $filterLine");

	system($filterLine);
	return $filteredDeltaFile;
}

sub _createMummerLine{
	my $self=shift;

	my $mummerLine = $self->mummerDirectory . 'nucmer --maxmatch';

	$mummerLine .= ' -b ' . $self->b if defined $self->b;
	$mummerLine .= ' -c ' . $self->c if defined $self->c;
	$mummerLine .= ' -d ' . $self->d if defined $self->d;
	$mummerLine .= ' -g ' . $self->g if defined $self->g;
	$mummerLine .= ' -l ' . $self->l if defined $self->l;
	$mummerLine .= ' -b ' . $self->b if defined $self->b;
	$mummerLine .= ' -p ' . $self->p if defined $self->p;
	$mummerLine .= ' ' . $self->referenceFile;
	$mummerLine .= ' ' . $self->queryFile;

	return $mummerLine;
}


sub _showCoords{
	my $self=shift;
	my $deltaFile = shift;

	#check for requirements
	unless(defined $self->coordsFile){
		$self->logger->logconfess("Coords file name required in _showCoords");
	}

	my $coordsLine = $self->mummerDirectory . 'show-coords -l -q -T ' . $deltaFile . ' > ' . $self->coordsFile;
	
	$self->logger->debug("Launching show-coords with $coordsLine");

	system($coordsLine);
}

1;
