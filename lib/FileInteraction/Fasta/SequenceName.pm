#!/usr/bin/perl

package FileInteraction::Fasta::SequenceName;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";


sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_initialize(@_);
    return $self;
}

sub arrayOfHeaders{
	my $self=shift;
	$self->{'_arrayOfHeaders'}=shift // return $self->{'_arrayOfHeaders'};
}

sub name{
	my $self=shift;
	$self->{'_name'}=shift // return $self->{'_name'};
}


##methods
sub _initialize{
	my $self=shift;

	if(@_){
		$self->getName(@_);
	}
}


sub getName{
	my($self)=shift;
	
	if(@_){
		my $name=shift;
		
		if($name =~ m/name=\|(\w+)\|/){
			$self->name($1);
		}
		elsif($name =~ m/lcl\|(\w+)\|/){
			$self->name($1);
		}
		elsif($name =~ m/(ref\|\w\w_\w\w\w\w\w\w|gb\|\w\w\w\w\w\w\w\w|emb\|\w\w\w\w\w\w\w\w)/){
			$self->name($1);
		}
		elsif($name =~ m/^(.+)\|Segment=/){
			$self->name($1);
		}
		elsif($name =~ m/^(.+)\|Length=/){
			$self->name($1);
		}
		else{
			$self->name($name);
		}
		$self->addToArrayOfHeaders($name);
	}
	else{
		print STDERR "no array sent!\n";
		exit(1);
	}	
}

sub addToArrayOfHeaders{
	my($self)=shift;
	
	if(@_){
		my $value=shift;
		my @array;
		@array = @{$self->arrayOfHeaders} if defined $self->arrayOfHeaders;
		push @array, $value;
		$self->arrayOfHeaders(\@array);
	}
	else{
		print STDERR 'nothing sent to addToArrayOfHeaders!' . "\n";
		exit(1);
	}
}

1;

