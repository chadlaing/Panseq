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
    }
    else{
    	#launch the panseq program
    	close STDERR;
    	close STDOUT;

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
	my $runMode = $q->param('runMode');

	my @loopVals;
	foreach my $val(@vals){
		unless ($self->_keepCurrentSetting($val, $runMode)){
			next;
		}

		# print STDERR "val: $val\nparamVal: " . $q->param($val) . "\nrunmode:$runMode\n";
		my %settings=(
			'VAL'=>$val // 0,
			'INFO'=>$q->param($val) // 0
		);
		push @loopVals, \%settings;
	}
	my %runMode=(
		'VAL'=>'runMode',
		'INFO'=>$runMode
	);
	push @loopVals, \%runMode;
	return \@loopVals;
}

sub _keepCurrentSetting{
	my $self=shift;
	my $key=shift;
	my $mode=shift;

	my %modeHash=(
		'email'=>1,
		'runMode'=>1
	);

	if($mode eq 'pan' || $mode eq 'novel'){
		$modeHash{'querySelected'}=1;
		$modeHash{'nucB'}=1;
		$modeHash{'nucC'}=1;
		$modeHash{'nucD'}=1;
		$modeHash{'nucG'}=1;
		$modeHash{'nucL'}=1;

		if($mode eq 'pan'){
			$modeHash{'coreGenomeThreshold'}=1;
			$modeHash{'percentIdentityCutoff'}=1;
			$modeHash{'fragmentationSize'}=1;
		}
		else{
			$modeHash{'referenceSelected'}=1;
		}
	}
	elsif($mode eq 'loci'){
		$modeHash{'numberOfLoci'}=1;
		$modeHash{'lociFile'}=1;
	}
	else{
		print STDERR "Unknown Panseq mode: $mode";
		exit(1);
	}


	if($modeHash{$key}){
		return 1;
	}
	else{
		return 0;
	}
}
1;




