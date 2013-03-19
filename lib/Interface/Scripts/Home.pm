#!/usr/bin/env perl
package Interface::Scripts::Home;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use parent 'Interface::Scripts::Panseq_Super';

sub setup{
	my $self=shift;


	$self->start_mode('display');
	$self->run_modes(
		'display'=>'display',
		'submit'=>'submit'
	);
}

sub display{
	my $self=shift;

	my $template = $self->load_tmpl('home.tmpl', die_on_bad_params=>0);
	return $template->output();
}

sub submit{
	my $self=shift;

	my $LOOP_VAL_REF = $self->_getFormSettings();

	my $template = $self->load_tmpl('submit.tmpl', die_on_bad_params=>0);
	$template->param(LOOP_VAL => $LOOP_VAL_REF);
	return $template->output();
}


sub _getFormSettings{
	my $self =shift;

	my $q = $self->query();
	my @vals = $q->param();

	my @loopVals;
	foreach my $val(@vals){
		my %settings=(
			'VAL'=>$val,
			'INFO'=>$q->param($val)
		);
		push @loopVals, \%settings;
	}
	return \@loopVals;
}


1;




