#!/usr/bin/env perl

package Roles::FlexiblePrinter;
#allows a unified, flexible printing strategy
use strict;
use warnings;
use FindBin;
use IO::File;
use lib "$FindBin::Bin/../";
use Role::Tiny;
 
sub outputFH{
	my $self =shift;
	
	if(@_){
		$self->{'_outputFilehandle'}=shift;
	}
	else{
		if(defined $self->{'_outputFilehandle'}){
			return $self->{'_outputFilehandle'};
		}
		elsif(defined $self->outputFile){
			$self->{'_outputFilehandle'}= (IO::File->new('>' . $self->outputFile) or die "Could not create FH in Roles::FlexiblePrinter");
			return $self->{'_outputFilehandle'};
		}
		else{
			return *STDOUT{IO};
		}
	}
}

#methods
sub printOut{
	my($self)=shift;
	
	my $outFH = $self->outputFH;
	print $outFH @_;
}

sub outputFile{
	my $self=shift;
	$self->{'_outputFile'}=shift // return $self->{'_outputFile'};
}

sub setNewOutputFile{
	my $self=shift;
	my $newFile=shift;

	if(defined $newFile){
		if($self->outputFH){
			$self->outputFH->close();
		}
		
		$self->outputFile($newFile);
		$self->outputFH(IO::File->new('>' . $newFile) or die "Could not create FH from $newFile in Roles::FlexiblePrinter");
	}
}

sub DESTROY{
	my $self=shift;

	if($self->outputFH){
		$self->outputFH->close();
	}
}

1;
