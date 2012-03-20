#!/usr/bin/perl
package BlastResultObject;

use FindBin::libs;
use MSA::BlastBased::BlastHitObject;
use FileInteraction::Fasta::SequenceName;

use Object::Tiny::RW qw{
  name
  firstHitName
  iteration_num
  query_def
  query_len
  hitHash
};

#methods

sub addHit {
	my ($self) = shift;

	if (@_) {
		my $hitName = shift;

		unless ( defined $self->hitHash && defined $self->hitHash->{$hitName} )
		{
			my $hit  = BlastHitObject->new();
			my $name = SequenceName->new($hitName);
			$self->firstHitName( $name->name )
			  unless defined $self->firstHitName
			;    #used when any hit object is needed, not a specific one
			$self->addToHash( 'hitHash', $name->name, $hit );
		}
	}
	else {
		print STDERR "no hit sent!\n";
		exit(1);
	}
}

sub addToHash {
	my ($self) = shift;

	if ( scalar(@_) >= 2 ) {
		my $hashName  = shift;
		my $hashValue = pop;

		my $hashRef = $self->$hashName;

		for ( my $i = 0 ; $i < scalar(@_) ; $i++ ) {
			my $key = $_[$i];

			if ( ( $i + 1 ) == scalar(@_) ) {
				$hashRef->{$key} = $hashValue;
			}
			else {
				$hashRef->{$key} = {} unless exists( $hashRef->{$key} );
				$hashRef = $hashRef->{$key};
			}
		}
		$self->$hashName($hashRef);

	}
	else {
		print STDERR "not enough values sent!";
		exit(1);
	}
}

1;
