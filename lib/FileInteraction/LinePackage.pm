#!/usr/bin/perl

package LinePackage;
use IO::File;

#object creation
use Object::Tiny::RW qw{
	line
	lineArray
};

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->initialize(@_);
    return $self;
}


##methods

sub stripLine{
	my($self)=shift;
	
	if(@_){
		my $line=shift;
		$line =~ s/[\n\f\r]//g;
		return $line;
	}
	else{
		print STDERR "nothing sent to be cleaned!\n";
		exit(1);
	}
}

sub initialize{
	my($self)=shift;
	
	if(@_){    	
    	my $delim;
    	if(scalar(@_)==2){
    		$delim=shift;
    	}
    	else{
    		$delim='\t';
    	}
    	$self->line($self->stripLine(@_));
    	$self->_createLineArray($delim);
    }
}

sub _createLineArray{
	my($self)=shift;
	
	if(@_){
		my $seperator=shift;
		
		my @lineArray = split($seperator,$self->line);
		$self->lineArray(\@lineArray);
	}
	else{
		print STDERR "no seperator specified!\n";
		exit(1);
	}
}

1;


