#!/usr/bin/perl

=pod

=head1 NAME

Example::Module::ClassTemplate - A Class that provides the following functionality:

=head1 SYNOPSIS


	use Example::Module::ClassTemplate;
	
	my $obj = Example::Module::ClassTemplate->new(
		'param1'=>'setting1',
		'param2'=>'setting2'
	);
	$obj->method1;
	$obj->method2;

=head1 DESCRIPTION

This module is explained here.

=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: 

=head1 AUTHOR

Your name (yourname@email.com)

=head2 Methods

=cut

package Modules::Fasta::MultiFastaSequenceName;

#includes
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Modules::Fasta::SequenceName;
use Log::Log4perl;

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


=head3 _initialize

Initializes the logger.
Assigns all values to class variables.
Anything else that the _initialize function does.

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 

    $self->logger->debug("Logger initialized in Modules::Fasta::MultiFastaSequenceName");  

    my %params = @_;

    #on object construction set all parameters
    foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			#logconfess calls the confess of Carp package, as well as logging to Log4perl
			$self->logger->logconfess("$key is not a valid parameter in Modules::NovelRegion::NovelRegionFinder");
		}
	}

	#required variables
	unless(defined $self->fileName){
		$self->logger->logconfess("fileName required in Modules::Fasta::MultiFastaSequenceName");
	}

	#process the file on object creation
	$self->_processMultiFastaFile();
}


=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}


=head3 fileName

The single file that contains multiple fasta sequences.

=cut

sub fileName{
	my $self=shift;
	$self->{'_fileName'}=shift // return $self->{'_fileName'};
}

=head3 sequenceNameHash

A hash that stores every fasta sequence in ->fileName in the following way:
	_sequenceNameHash->{Modlues::Fasta::SequenceName->name} = Modlues::Fasta::SequenceName
This allows access to all the fastaHeaders for a given Modlues::Fasta::SequenceName via the
Modlues::Fasta::SequenceName->arrayOfHeaders method

=cut

sub genomeNameFromContig{
	my $self=shift;
	$self->{'__sequenceNameHash'}=shift // return $self->{'__sequenceNameHash'};
}

sub contigNamesFromGenome{
	my $self=shift;
	$self->{'_contigNamesFromGenome'}=shift // return $self->{'_contigNamesFromGenome'};
}

sub orderedGenomeNames{
	my $self=shift;
	$self->{'_orderedGenomeNames'}=shift // return $self->{'_orderedGenomeNames'};
}

sub _processMultiFastaFile{
	my $self=shift;

	my %genomeNames;
	my %contigToGenome;
	my %genomeToContigs;
	my $inFH = IO::File->new('<' . $self->fileName) or die "$!";
	while(my $line = $inFH->getline()){
		unless($line =~ m/^>/){
			next;
		}
		$line =~ s/\R//g;
		$line =~ s/>//;
		my $sn = Modules::Fasta::SequenceName->new($line);
		$genomeNames{$sn->name}=1;
		$contigToGenome{$line}=$sn->name;
		push @{$genomeToContigs{$sn->name}}, $line;
	}
	$inFH->close();

	$self->orderedGenomeNames([sort keys %genomeNames]);
	$self->genomeNameFromContig(\%contigToGenome);
	$self->contigNamesFromGenome(\%genomeToContigs);
}


1;


