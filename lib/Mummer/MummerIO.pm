#!/usr/bin/perl

package MummerIO;

use strict;
use warnings;
use IO::File;
use POSIX;

#object creation
use Object::Tiny::RW qw{
	mummerDirectory
	deltaFile
	b
	c
	d
	g
	l
	p
	nooptimize
	nosimplify
	noextend
};

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_mummerIOInitialize(@_);
	return $self;
}

##methods
sub _mummerIOInitialize{
	my($self)=shift;
	
	if(@_){
		my $programDir = shift;
		$self->mummerDirectory($programDir);
	}
	else{
		print STDERR "MUMmer requires the program directory to execute!\n";
		exit(1);
	}
}

sub runMummer {
	my($self) = shift;
	
	if(@_){
		my %params = @_;
		
		my $systemLine = $self->mummerDirectory . 'nucmer --maxmatch';		
		
		#check for previously set values
		if(defined $self->b){$params{'b'}=$self->b}
		if(defined $self->c){$params{'c'}=$self->c}
		if(defined $self->d){$params{'d'}=$self->d}
		if(defined $self->g){$params{'g'}=$self->g}
		if(defined $self->l){$params{'l'}=$self->l}
		if(defined $self->p){$params{'p'}=$self->p}
		if(defined $self->nooptimize){$params{'nooptimize'}=1}
		if(defined $self->nosimplify){$params{'nosimplify'}=1}
		if(defined $self->noextend){$params{'noextend'}=1}		
		
		#check/override values with specific run values
		if (defined $params{'b'}) { $systemLine .= ' -b ' . $params{'b'}}
		if (defined $params{'c'}) { $systemLine .= ' -c ' . $params{'c'}}
		if (defined $params{'d'}) { $systemLine .= ' -d ' . $params{'d'}}
		if (defined $params{'g'}) { $systemLine .= ' -g ' . $params{'g'}}
		if (defined $params{'l'}) { $systemLine .= ' -l ' . $params{'l'}}
		if (defined $params{'p'}) { $systemLine .= ' -p ' . $params{'p'}}
		if (defined $params{'nooptimize'}) { $systemLine .= ' --nooptimize'}
		if (defined $params{'nosimplify'}) { $systemLine .= ' --nosimplify'}
		if (defined $params{'noextend'}) { $systemLine .= ' --noextend'}
		
		if (defined $params{'referenceFile'}) { $systemLine .= ' ' . $params{'referenceFile'}}
		if (defined $params{'queryFile'}) { $systemLine .= ' ' . $params{'queryFile'}}
		
		#store delta file
		my $prefix = $params{'p'} || 'out';
		$self->deltaFile($prefix . 'delta');

		system($systemLine);
	}
	else {
		print STDERR
		  "no query and/or reference file or runDirectory defined!\n";
		exit(1);
	}
}

1;
