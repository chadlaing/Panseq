#!/usr/bin/perl

package SequenceName;

#object creation
use Object::Tiny::RW qw{
	name
	arrayOfHeaders
};


sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->getName(@_) if @_;
    return $self;
}

##methods
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

