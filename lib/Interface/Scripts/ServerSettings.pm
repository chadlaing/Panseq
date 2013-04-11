#!/usr/bin/env perl

package Interface::Scripts::ServerSettings;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";

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

    $self->_getSettingsFromConfigurationFile(@_);
}

#class variables
sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub emailDefaultsFile{
	my $self=shift;
	$self->{'_emailDefaults'} = shift // return $self->{'_emailDefaults'};
}

sub baseDirectory{
	my $self=shift;
	$self->{'_baseDirectory'}=shift // return $self->{'_baseDirectory'};
}

sub queryDirectory{
	my $self=shift;
	$self->{'_queryDirectory'}=shift // return $self->{'_queryDirectory'};
}

sub referenceDirectory{
	my $self=shift;
	$self->{'_referenceDirectory'}=shift // return $self->{'_referenceDirectory'};
}

sub blastDirectory{
	my $self=shift;
	$self->{'_blastDirectory'}=shift // return $self->{'_blastDirectory'};
}

sub mummerDirectory{
	my $self=shift;
	$self->{'_mummerExecutable'}=shift // return $self->{'_mummerExecutable'};
}

sub muscleExecutable{
	my $self=shift;
	$self->{'_muscleExecutable'}=shift // return $self->{'_muscleExecutable'};
}

sub fastaFileDirectory{
	my $self=shift;
	$self->{'_fastaFileDirectory'}=shift // return $self->{'_fastaFileDirectory'};
}

sub numberOfCores{
	my $self=shift;
	$self->{'_numberOfCores'}=shift // return $self->{'_numberOfCores'};
}

sub outputDirectory{
	my $self=shift;
	$self->{'_outputDirectory'}=shift // return $self->{'_outputDirectory'};
}


sub _getSettingsFromConfigurationFile{
	my $self=shift;

	if (@_) {
		my $fileName = shift;
		my $inFile = IO::File->new( '<' . $fileName ) or die "\nCannot open $fileName\n $!";
		
		while ( my $line = $inFile->getline ) {
			next unless $line =~ m/\t/;
			
			$line =~ s/\R//g;
			my @la = split(/\t/,$line);

			my $setting = $la[0] // undef;
			my $value   = $la[1] // undef;
			
			if($self->can($setting)){
				$self->$setting($value);
			}
			else{
				$self->logger->fatal("Unknown setting $setting");
				exit(1);
			}
		}
		$inFile->close();
	}
	else {
		$self->logger->fatal("No configuration file specified");
		exit(1);
	}
}


1;
