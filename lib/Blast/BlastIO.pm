#!/usr/bin/perl

package Blast::BlastIO;
use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin";

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


sub programDirectory{
	my $self=shift;
	$self->{'_programDirectory'}=shift // return $self->{'_programDirectory'};
}

sub type{
	my $self=shift;
	$self->{'_type'}=shift // return $self->{'_type'};
}

sub task{
	my $self=shift;
	$self->{'_task'}=shift // return $self->{'_task'};
}

sub dust{
	my $self=shift;
	$self->{'_dust'}=shift // return $self->{'_dust'};
}

sub query{
	my $self=shift;
	$self->{'_query'}=shift // return $self->{'_query'};
}

sub db{
	my $self=shift;
	$self->{'_db'}=shift // return $self->{'_db'};
}

sub out{
	my $self=shift;
	$self->{'_out'}=shift // return $self->{'_out'};
}

sub evalue{
	my $self=shift;
	$self->{'_evalue'}=shift // return $self->{'_evalue'};
}

sub word_size{
	my $self=shift;
	$self->{'_word_size'}=shift // return $self->{'_word_size'};
}

sub outfmt{
	my $self=shift;
	$self->{'_outfmt'}=shift // return $self->{'_outfmt'};
}

sub num_descriptions{
	my $self=shift;
	$self->{'_num_descriptions'}=shift // return $self->{'_num_descriptions'};
}

sub num_alignments{
	my $self=shift;
	$self->{'_num_alignments'}=shift // return $self->{'_num_alignments'};
}

sub max_target_seqs{
	my $self=shift;
	$self->{'_max_target_seqs'}=shift // return $self->{'_max_target_seqs'};
}

sub num_threads{
	my $self=shift;
	$self->{'_num_threads'}=shift // return $self->{'_num_threads'};
}

sub best_hit_score_edge{
	my $self=shift;
	$self->{'_best_hit_score_edge'}=shift // return $self->{'_best_hit_score_edge'};
}

sub best_hit_overhang{
	my $self=shift;
	$self->{'_best_hit_overhang'}=shift // return $self->{'_best_hit_overhang'};
}

#methods
sub _initialize{
	my($self)=shift;
	
	my $paramsRef=shift;	
	$self->setValuesFromHash($paramsRef);	
}


sub runBlastn{
	my($self)=shift;
	
	if(1){
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
		$systemLine .= ' -best_hit_score_edge ' . $self->best_hit_score_edge if (defined $self->best_hit_score_edge);
		$systemLine .= ' -best_hit_overhang ' . $self->best_hit_overhang if (defined $self->best_hit_overhang);

		
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
	
	if(1){
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
		$self->best_hit_score_edge($hashRef->{'best_hit_score_edge'}) if defined $hashRef->{'best_hit_score_edge'};
		$self->best_hit_overhang($hashRef->{'best_hit_overhang'}) if defined $hashRef->{'best_hit_overhang'};
		
		return;
	}
	else{
		print STDERR "nothing sent to setValuesFromHash\n";
		exit(1);
	}
}


1;