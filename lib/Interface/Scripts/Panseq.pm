#!/usr/bin/env perl
package Interface::Scripts::Panseq;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use parent 'Interface::Scripts::Panseq_Super';

sub setup{
	my $self=shift;


	$self->start_mode('display');
	$self->run_modes(
		'analyses'=>'analyses',
		'submit'=>'submit',
		'home'=>'home',
		'contact'=>'contact'
	);

	$self->batchFile('/home/chad/batfile.txt');
}

sub batchFile{
	my $self=shift;
	$self->{'_batchFile'}=shift // return $self->{'_batchFile'};
}

sub home{
	my $self=shift;

	my $template=$self->load_tmpl('home.tmpl',die_on_bad_params=>0);
	return $template->output();
}

sub analyses{
	my $self=shift;

	my $template = $self->load_tmpl('analyses.tmpl', die_on_bad_params=>0);
	return $template->output();
}

sub contact{
	my $self=shift;

	my $template = $self->load_tmpl('contact.tmpl', die_on_bad_params=>0);
	return $template->output();
}

sub submit{
	my $self=shift;

	my $pid = fork();
    if(!defined $pid){
            die "cannot fork process!\n $!";
    };


    if($pid){
        #display alldone!
        my $LOOP_VAL_REF = $self->_getFormSettings();
		my $template = $self->load_tmpl('submit.tmpl', die_on_bad_params=>0);
		$template->param(LOOP_VAL => $LOOP_VAL_REF);
		return $template->output();
		exit(0);
    }
    else{
    	#launch the panseq program
    	$self->_launchPanseq();
	}

	
}


sub _launchPanseq{
	my $self=shift;

	$self->_createBatchFile();
}


sub _createBatchFile{
	my $self=shift;

	my $q = $self->query();
	my @params = $q->param();

	my $batchFH = IO::File->new('>' . $self->batchFile) or die "$!";

	foreach my $setting(@params){
		$batchFH->print($setting . "\t" . $q->param($setting) . "\n");
	}

	$batchFH->close();
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




