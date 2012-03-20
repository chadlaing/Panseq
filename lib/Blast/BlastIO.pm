#!/usr/bin/perl

package BlastIO;
use strict;
use warnings;
use Carp;

#commenting for the comment lovers
use Object::Tiny::RW qw{
	programDirectory
	type
	task
	dust
	query
	db
	out
	evalue
	word_size
	outfmt
	num_descriptions
	num_alignments
	max_target_seqs
	num_threads	
};

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_blastIOInitialize(@_);
	return $self;
}

#methods
sub _blastIOInitialize{
	my($self)=shift;
	
	my $paramsRef=shift;	
	$self->setValuesFromHash($paramsRef);	
}


sub runBlastn{
	my($self)=shift;
	
	if(@_){
		my %settingsHash=@_;
		
		$self->setValuesFromHash(\%settingsHash);
		
		unless((defined $self->db) && (defined $self->query)){
			print STDERR "runBlastn requires at minimum db and query to be defined!\n";
			exit(1);
		}
		$self->programDirectory // confess ('programDirectory is undefined');
		$self->type // confess ('type is undefined');
		
		my $systemLine=$self->programDirectory . $self->type;
		
		$systemLine .= ' -task ' . $self->task if (defined $self->task);
		$systemLine .= ' -dust ' . $self->dust if (defined $self->dust);
		$systemLine .= ' -query ' . $self->query if (defined $self->query);
		$systemLine .= ' -db ' . $self->db if (defined $self->db);
		$systemLine .= ' -out ' . $self->out if (defined $self->out);
		$systemLine .= ' -evalue ' . $self->evalue if (defined $self->evalue);
		$systemLine .= ' -word_size ' . $self->word_size if (defined $self->word_size);
		$systemLine .= ' -outfmt ' . $self->outfmt if (defined $self->outfmt);
		$systemLine .= ' -num_descriptions ' . $self->num_descriptions if (defined $self->num_descriptions);
		$systemLine .= ' -num_alignments ' . $self->num_alignments if (defined $self->num_alignments);
		$systemLine .= ' -max_target_seqs ' . $self->max_target_seqs if (defined $self->max_target_seqs);
		$systemLine .= ' -num_threads ' . $self->num_threads if (defined $self->num_threads);
		
		system($systemLine);
		return;
	}
	else{
		print STDERR "incorrect number of arguments sent to runBlastn\n";
		exit(1);
	}
}

sub setValuesFromHash{
	my($self)=shift;
	
	if(@_){
		my $hashRef=shift;
		
		$self->programDirectory($hashRef->{'blastDirectory'}) if defined $hashRef->{'blastDirectory'};
		$self->type($hashRef->{'type'}) if defined $hashRef->{'type'} ;
		$self->task($hashRef->{'task'}) if defined $hashRef->{'task'} ;
		$self->dust($hashRef->{'dust'}) if defined $hashRef->{'dust'};
		$self->query($hashRef->{'query'}) if defined $hashRef->{'query'};
		$self->db($hashRef->{'db'}) if defined $hashRef->{'db'};
		$self->out($hashRef->{'out'}) if defined $hashRef->{'out'};
		$self->evalue($hashRef->{'evalue'}) if defined $hashRef->{'evalue'};
		$self->word_size($hashRef->{'word_size'}) if defined $hashRef->{'word_size'};
		$self->outfmt($hashRef->{'outfmt'}) if defined $hashRef->{'outfmt'};
		$self->num_descriptions($hashRef->{'num_descriptions'}) if defined $hashRef->{'num_descriptions'};
		$self->num_alignments($hashRef->{'num_alignments'}) if defined $hashRef->{'num_alignments'};
		$self->max_target_seqs($hashRef->{'max_target_seqs'}) if defined $hashRef->{'max_target_seqs'};
		$self->num_threads($hashRef->{'num_threads'}) if defined $hashRef->{'num_threads'};
		
		return;
	}
	else{
		print STDERR "nothing sent to setValuesFromHash\n";
		exit(1);
	}
}


1;