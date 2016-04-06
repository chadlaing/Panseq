#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use CGI;

my $cgi = CGI->new();
my $pid = fork();
if(!defined $pid){
    die "cannot fork process!\n $!";
};

if($pid){
    #the redirect url goes here
    print $cgi->redirect("/page/index.html");

}
else{
    print STDERR "processing the submission\n";
#    #launch the panseq program
#    close STDERR;
#    close STDOUT;
#    my $q=$self->query();
#
#    #    	if($q->param('runMode') eq 'loci'){
#    #    		$self->_launchLociSelector($q);
#    #    	}
#    #    	else{
#    $self->_launchPanseq();
#    #    	}
}