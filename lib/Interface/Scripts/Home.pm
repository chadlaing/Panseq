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
		'display'=>'display'
	);
}

sub display{
	my $self=shift;

	my $template = $self->load_tmpl('home.tmpl', die_on_bad_params=>0);
	return $template->output();
}

1;


