#!/usr/bin/env perl
package Modules::Fasta::FastaFileSplitter;


=pod

=head1 NAME

Modules::Fasta::FastaFileSplitter - Splits a single fasta file into a given number of ~equal sized output fasta files.

=head1 SYNOPSIS

Given a number of files to create, we iterate through the master file, adding one fasta sequence to each of
the output files in turn. The output is accumulated in an array, and then printed to the files at the end.
We add sequences in this iterative fashion to help in the downstream BLAST step, where "core" regions tend to
cluster, making some output files significantly larger than others. This leads to most cores sitting idle while 
a few work to finish processing the largest blast output. By adding sequence iteratively to the output files, the BLAST
output should be relatively homogeneous.

=head1 DESCRIPTION



=head2 Methods


=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing gmail com)

=cut

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Carp;
use IO::File;
use Log::Log4perl;


#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

#class variables
sub logger{
	my $self=shift;
	$self->{'_FastaFileSplitter_logger'} = shift // return $self->{'_FastaFileSplitter_logger'};
}

sub inputFile{
	my $self=shift;
	$self->{'_inputFile'}=shift // return $self->{'_inputFile'};
}

sub numberOfSplits{
	my $self=shift;
	$self->{'_numberOfSplits'}=shift // return $self->{'_numberOfSplits'};
}

sub baseDirectory{
    my $self = shift;
    $self->{'_baseDirectory'} = shift // return $self->{'_baseDirectory'};   
}


#methods
sub _initialize {
	my $self = shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Fasta::FastaFileSplitter");


	my %params = @_;
    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::Fasta::FastaFileSplitter");
		}
	}	

	unless(defined $self->inputFile){
		$self->logger->logconfess("Modules::Fasta::FastaFileSplitter requires:\n
			\t'inputFile'\n
			to be defined. One or more parameters not defined.
		");
	}

	unless(defined $self->numberOfSplits){
		$self->logger->info("Number of splits not specified, will not split input file.");
		$self->numberOfSplits(1);
	}

	unless(defined $self->baseDirectory){
		$self->logger->logconfess("baseDirectory is not specified");
	}
}

sub splitFastaFile {
	my $self = shift;
	
	my @splitFiles;
	$self->logger->info("Splitting " . $self->inputFile . " into " . $self->numberOfSplits . " files");

	my @splitArray;

	if($self->numberOfSplits == 1){
		$splitFiles[0] = $self->inputFile;
	}
	else{
		
		my $inFH = IO::File->new('<' . $self->inputFile) or die "$!";
		my $outputPrefix = 0;

		while(my $line = $inFH->getline()){
			
			if($line =~ m/^>/){
				$outputPrefix++;
				if($outputPrefix > $self->numberOfSplits){
					$outputPrefix=1;
				}					
			}
			if(defined $splitArray[$outputPrefix]){
				$splitArray[$outputPrefix] .= $line;
			}
			else{
				$splitArray[$outputPrefix] = $line;
			}
		}		
		$inFH->close();
		#print out the array values to the separate files
		for my $i(1..$self->numberOfSplits){
			my $outFileName = $self->baseDirectory . $i . '_split.fasta';
			my $outFH = IO::File->new('>' . $outFileName) or die "$!";
			$outFH->print($splitArray[$i]);
			push @splitFiles, $outFileName;
			$outFH->close();
		}
	}		
	return \@splitFiles;
}

1;
