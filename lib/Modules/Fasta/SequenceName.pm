#!/usr/bin/perl

=pod

=head1 NAME

Modules::Fasta::SequenceName - Takes a fasta header and returns the genome name. Stores all headers from a multi-fasta file.

=head1 SYNOPSIS


	use Modules::Fasta::SequenceName;
	
	my $obj = Modules::Fasta::SequenceName->new(
		'gi|148566106|gb|CP000711.1| Enterobacteria phage CUS-3, complete genome'
	);
	$obj->name;
	$obj->arrayOfHeaders;

=head1 DESCRIPTION

This module is used to unify the name of a genome when given a multi-fasta sequence, such as the case for a draft genome.
It allows the storing of multiple headers that all share the same "name".
Removes the '>' fasta header ID by default.

=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: 

=head1 AUTHOR

Your name (chadlaing@gmail.com)

=head2 Methods

=cut

package Modules::Fasta::SequenceName;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";

use Log::Log4perl;
use Carp; #this enables a stack trace unlike warn/die

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

    #init anonymous array
    #this has to be before it us used
    #so putting it at the end of the init sub is a bad idea
	$self->arrayOfHeaders([]);

    #on object construction set all parameters
	if(scalar(@_)==1){
		$self->name($self->_getName(@_));
	}
	else{
		#logconfess calls the confess of Carp package, as well as logging to Log4perl
		$self->logger->logconfess("Modules::Fasta::SequenceName constructor requires a single value");
	}

}

=head3 logger

Stores an instance of the logger.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}


=head3 arrayOfHeaders

All fasta headers sharing the same name can be stored here.
Although only the header specified in new() is set initially,
the arrayOfHeaders is a public method that could be used outside of the Class
to store all the headers of a multi-fasta file.

=cut

sub arrayOfHeaders{
	my $self=shift;
	$self->{'_arrayOfHeaders'}=shift // return $self->{'_arrayOfHeaders'};
}

=head3 name

The method that returns the genome "name" from the header passed to new.
Only the initialization of the object is allowed to set the name.
After that, only headers can be added to arrayOfHeaders, but the name is permanent.

=cut

sub name{
	my $self=shift;
	$self->{'_name'}=shift // return $self->{'_name'};
}


=head3 _getName

Private method that does the extracting of the "name" from the fasta header passed to new().
In order it looks for name=||, lcl||, ref||, gb||, emb||, dbj||, gi||, |Segment=, |Length=.
The name between the bars for name and lcl is used. For all other names except Segment= and Length=, the identifier
is included. Eg. the name returned is ref|NC_192832. For Segment and Length, the entire header up to the name is used.
If no match occurs for the name, the entire fasta header is returned.

=cut

sub _getName{
	my $self =shift;
	my $originalName=shift;
	
	my $newName;
	if($originalName =~ m/name=\|(\w+)\|/){
		$newName = $1;
	}
	elsif($originalName =~ m/lcl\|([\w-]*)\|/){
		$newName = $1;
	}
	elsif($originalName =~ m/(ref\|\w\w_\w\w\w\w\w\w|gb\|\w\w\w\w\w\w\w\w|emb\|\w\w\w\w\w\w\w\w|dbj\|\w\w\w\w\w\w\w\w)/){
		$newName = $1;
	}
	elsif($originalName =~ m/(gi\|\d+)\|/){
		$newName = $1;
	}
	elsif($originalName =~ m/^(.+)\|Segment=/){
		$newName = $1;
	}
	elsif($originalName =~ m/^(.+)\|Length=/){
		$newName = $1;
	}
	else{
		$newName = $originalName;
	}
	$newName =~ s/^>//;
	push @{$self->arrayOfHeaders},$originalName;
	return $newName;
}

1;

