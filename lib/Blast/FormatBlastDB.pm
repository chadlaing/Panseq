#!/usr/bin/perl

package FormatBlastDB;

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

#methods
sub initialize{
	my($self)=shift;
	
	if(@_){
		my $dir=shift;
		$self->programDirectory($dir);
	}
	else{
		print STDERR "Program directory must be specified in FormatBlastDB\n";
		exit(1);
	}
	
	#logging
	$self->logger( Log::Log4perl->get_logger());
}


sub runMakeBlastDb{
	my($self)=shift;
	
	if(@_){
		my %runValues=@_;	
		
		unless((defined $runValues{'in'})){
			print STDERR "Input file must be identified in runMakeBlastDb\n";
			exit(1);
		}
		
		my $systemLine =$self->programDirectory . 'makeblastdb -dbtype ' . $runValues{'dbtype'};
		
		$systemLine .= ' -out '. $runValues{'out'} if defined $runValues{'out'};
		$systemLine .= ' -title ' . $runValues{'title'} if defined $runValues{'title'};
		$systemLine .= ' -in ' . $runValues{'in'} if defined $runValues{'in'};
		$systemLine .= ' -logfile ' . $runValues{'logFile'} if defined $runValues{'logFile'};
		
		$self->logger->info("Running makeblastdb with the following: $systemLine");
		system($systemLine); 
	}
	else{
		print STDERR "formatdb missing parameter hash!\n";
		exit(1);
	}	
}

1;
