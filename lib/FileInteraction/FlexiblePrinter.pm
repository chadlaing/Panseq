#!/usr/bin/perl

package FlexiblePrinter;
#allows a unified, flexible printing strategy

sub new{
	my($class)  = shift;
    my $self= {};
    bless ($self, $class);
    $self->_flexiblePrinterInitialize(@_);  
    return $self;
}


#get/set
sub _flexiblePrinterInitialize{
	my $self=shift;
	
	if(@_){
    	my $handle=shift;
    	$self->outputFilehandle($handle);
    } 
}

 
sub outputFilehandle{
	my($self)=shift;
	
	if(@_){
		$self->{'_outputFilehandle'}=shift;
	}
	else{
		return $self->{'_outputFilehandle'} || *STDOUT{IO};
	}
}

#methods
sub printOut{
	my($self)=shift;
	
	my $outFile = $self->outputFilehandle;
	print $outFile @_;
}

1;
