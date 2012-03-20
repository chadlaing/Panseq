#!/usr/bin/perl

#10-5-11 - Script intended for use with blastn 2.2.
#@TODO:Merge with BlastIO 28-10-11
#functionality of blastn included:
#db
#db_soft_mask
#subject
#query
#out
#evalue
#word_size
#gapoen
#gapextend
#perc_identity
#penalty
#min_raw_gapped_score
#num_threads
#num_descriptions
#num_alingments
#outfmt
#max_target_seq

package BlastnCmd;

use warnings;
use strict;
use diagnostics;
use Carp;
use Scalar::Util qw(looks_like_number);

my ($isDbUsed)           = 0;
my ($isSubjectUsed)      = 0;
my ($isDb_soft_maskUsed) = 0;

#constructor
sub new {
	my ($class) = shift;
	my ($self)  = {};
	bless $self, $class;
	return $self;
}

#set functions

#file path directory of blastn
sub setDirectory {
	my ( $self, $directory ) = @_;
	$self->{_directory} = $directory if defined($directory);
}

#include file path
sub setQuery {
	my ( $self, $query ) = @_;
	$self->{_query} = $query if defined($query);
}

#include file path
sub setSubject {
	if ($isDbUsed) {
		confess
"Already started use of -db and cannot also use -subject in blastn call";
	}
	elsif ($isDb_soft_maskUsed) {
		confess
"Already started use of -db_soft_mask and cannot also use -subject in blastn call";
	}
	else {
		$isSubjectUsed = 1;
		my ( $self, $subject ) = @_;
		$self->{_subject} = $subject if defined($subject);
	}
}

#include file path

sub setDb {
	if ($isSubjectUsed) {
		confess
"Already started use of -subject and cannot also use -db in blastn call";
	}
	else {
		$isDbUsed = 1;
		my ( $self, $db ) = @_;
		$self->{_db} = $db if defined($db);
	}
}

sub setWord_size {
	my ( $self, $word_size ) = @_;
	if ( ( $word_size =~ m/^\d+$/ ) && ( $word_size >= 4 ) ) {
		$self->{_word_size} = $word_size if defined($word_size);
	}
	else {
		confess "-word_size value must and integer greater or equal to 4";
	}

}

#include file path and approproite file extension type
sub setOut {
	my ( $self, $out ) = @_;
	$self->{_out} = $out if defined($out);
}

sub setOutfmt {
	my ( $self, $outfmt ) = @_;
	if ( ( $outfmt =~ /^-?\d+$/ ) && ( $outfmt < 12 ) && ( $outfmt >= 0 ) ) {
		$self->{_outfmt} = $outfmt if defined($outfmt);
	}
	else {
		confess "-outfmt input values must be integers between 0 and 11";
	}
}

sub setDb_soft_mask {
	my ( $self, $db_soft_mask ) = @_;
	if ($isSubjectUsed) {
		confess
"Already started use of -subject and cannot also use -db_soft_mask in blastn call";
	}
	else {
		$isDb_soft_maskUsed = 1;
		$self->{_db_soft_mask} = $db_soft_mask if defined($db_soft_mask);
	}
}

sub setEvalue {
	my ( $self, $evalue ) = @_;
	if ( ( looks_like_number($evalue) && ( $evalue > 0 ) ) ) {
		$self->{_evalue} = $evalue if defined($evalue);
	}
	else {
		confess "-evalue threshold must be an number greater than 0";
	}
}

sub setGapopen {
	my ( $self, $gapopen ) = @_;
	if ( ( $gapopen =~ /^-?\d+$/ ) && ( $gapopen >= 0 ) ) {
		$self->{_gapopen} = $gapopen if defined($gapopen);
	}
	else {
		confess
		  "-gapopen penalty value must be an integer greater or equal to zero";
	}
}

sub setGapextend {
	my ( $self, $gapextend ) = @_;
	if ( ( $gapextend =~ /^-?\d+$/ ) && ( $gapextend > 0 ) ) {
		$self->{_gapextend} = $gapextend if defined($gapextend);
	}
	else {
		confess "-gapextend penalty value must be an integer greater than zero";
	}
}

sub setPerc_identity {
	my ( $self, $perc_identity ) = @_;
	if (   looks_like_number($perc_identity)
		&& ( $perc_identity >= 0 )
		&& ( $perc_identity <= 100 ) )
	{
		$self->{_perc_identity} = $perc_identity if defined($perc_identity);
	}
	else {
		confess "-perc_identity a number must greater than zero";
	}

}

sub setPenalty {
	my ( $self, $penalty ) = @_;
	if ( ( $penalty =~ /^-?\d+$/ ) && ( $penalty <= 0 ) ) {
		$self->{_penalty} = $penalty if defined($penalty);
	}
	else {
		confess
		  "-penalty base mismatch value must be an integer less than zero";
	}
}

sub setMin_raw_gapped_score {
	my ( $self, $min_raw_gapped_score ) = @_;
	if ( $min_raw_gapped_score =~ /^-?\d+$/ )

	{
		$self->{_min_raw_gapped_score} = $min_raw_gapped_score
		  if defined($min_raw_gapped_score);
	}
	else {
		confess "-min_raw_gapped_score threshold must be an integer";
	}
}

sub setNum_threads {
	my ( $self, $num_threads ) = @_;
	if ( ( $num_threads =~ /^-?\d+$/ ) && ( $num_threads > 0 ) ) {
		$self->{_num_threads} = $num_threads if defined($num_threads);
	}
	else {
		confess "-num_threads value must be an integer greater than zero";
	}
}

sub setNum_descriptions {
	my ( $self, $num_descriptions ) = @_;
	if ( ( $num_descriptions =~ /^-?\d+$/ ) && ( $num_descriptions >= 0 ) ) {
		$self->{_num_descriptions} = $num_descriptions
		  if defined($num_descriptions);
	}
	else {
		confess
		  "-num_descriptions value must be an integer greater or equal to zero";
	}
}

sub setNum_alignments {
	my ( $self, $num_alignments ) = @_;
	if ( ( $num_alignments =~ /^-?\d+$/ ) && ( $num_alignments >= 0 ) ) {
		$self->{_num_alignments} = $num_alignments
		  if defined($num_alignments);
	}
	else {
		confess
		  "-num_alignments value must be an integer greater or equal to zero";
	}
}

sub setMax_target_seqs {
	my ( $self, $max_target_seqs ) = @_;
	if ( ( $max_target_seqs =~ /^-?\d+$/ ) && ( $max_target_seqs > 0 ) ) {
		$self->{_max_target_seqs} = $max_target_seqs
		  if defined($max_target_seqs);
	}
	else {
		confess "-max_target_seqs value must be an integer greater zero";
	}
}

#get functions

sub getDirectory {
	my ($self) = @_;
	return $self->{_directory};
}

sub getQuery {
	my ($self) = @_;
	return $self->{_query};
}

sub getSubject {
	my ($self) = @_;
	return $self->{_subject};
}

sub getDb {
	my ($self) = @_;
	return $self->{_db};
}

sub getWord_size {
	my ($self) = @_;
	return $self->{_word_size};
}

sub getOut {
	my ($self) = @_;
	1;
	return $self->{_out};
}

sub getOutfmt {
	my ($self) = @_;
	return $self->{_outfmt};
}

sub getDb_soft_mask {
	my ($self) = @_;
	return $self->{_db_soft_mask};
}

sub getEvalue {
	my ($self) = @_;
	return $self->{_evalue};
}

sub getGapopen {
	my ($self) = @_;
	return $self->{_gapopen};
}

sub getGapextend {
	my ($self) = @_;
	return $self->{_gapextend};
}

sub getPerc_identity {
	my ($self) = @_;
	return $self->{_perc_identity};
}

sub getPenalty {
	my ($self) = @_;
	return $self->{_penalty};
}

sub getMin_raw_gapped_score {
	my ($self) = @_;
	return $self->{_min_raw_gapped_score};
}

sub getNum_threads {
	my ($self) = @_;
	return $self->{_num_threads};
}

sub getNum_descriptions {
	my ($self) = @_;
	return $self->{_num_descriptions};
}

sub getNum_alignments {
	my ($self) = @_;
	return $self->{_num_alignments};
}

sub getMax_target_seqs {
	my ($self) = @_;
	return $self->{_max_target_seqs};
}

#primary blast method

sub runBlastn {
	my ($self) = shift;
	if ( $self->getQuery() ) {
		my ($systemLine) =
		  $self->getDirectory() . 'blastn -query ' . $self->getQuery();
		if ( $isDbUsed || $isSubjectUsed ) {

			if ( $self->getDb() ) {
				$systemLine = $systemLine . ' -db ' . $self->getDb();
			}
			if ( $self->getSubject() ) {
				$systemLine = $systemLine . ' -subject ' . $self->getSubject();
			}
			if ( $self->getWord_size() ) {
				$systemLine =
				  $systemLine . ' -word_size ' . $self->getWord_size();
			}
			if ( $self->getOut() ) {
				$systemLine = $systemLine . ' -out ' . $self->getOut();
			}
			if ( $self->getOutfmt() ) {
				$systemLine = $systemLine . ' -outfmt ' . $self->getOutfmt();
			}
			if ( $self->getDb_soft_mask() ) {
				$systemLine =
				  $systemLine . ' -db_soft_mask ' . $self->getDb_soft_mask();
			}
			if ( $self->getEvalue() ) {
				$systemLine = $systemLine . ' -evalue ' . $self->getEvalue();
			}
			if ( $self->getGapopen() ) {
				$systemLine = $systemLine . ' -gapopen ' . $self->getGapopen();
			}
			if ( $self->getGapextend() ) {
				$systemLine =
				  $systemLine . ' -gapextend ' . $self->getGapextend();
			}
			if ( $self->getPerc_identity() ) {
				$systemLine =
				  $systemLine . ' -perc_identity ' . $self->getPerc_identity();
			}
			if ( $self->getPenalty() ) {
				$systemLine = $systemLine . ' -penalty ' . $self->getPenalty();
			}
			if ( $self->getMin_raw_gapped_score() ) {
				$systemLine =
				    $systemLine
				  . ' -min_raw_gapped_score '
				  . $self->getMin_raw_gapped_score();
			}
			if ( $self->getNum_threads() ) {
				$systemLine =
				  $systemLine . ' -num_threads ' . $self->getNum_threads();
			}
			if ( $self->getNum_descriptions() ) {
				$systemLine =
				    $systemLine
				  . ' -num_descriptions '
				  . $self->getNum_descriptions();
			}
			if ( $self->getNum_alignments() ) {
				$systemLine =
				    $systemLine
				  . ' -num_alignments '
				  . $self->getNum_alignments();
			}
			if ( $self->getMax_target_seqs() ) {
				$systemLine =
				    $systemLine
				  . ' -max_target_seqs '
				  . $self->getMax_target_seqs();
			}
			system($systemLine);
		}
		else { confess "no -database or -subject value used" }
	}
	else {
		confess "No -query value used";
	}
}
1;
