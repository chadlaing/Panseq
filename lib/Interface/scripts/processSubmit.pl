#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use CGI;
use Data::Dumper;

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

    #we need to determine what run mode (pan, novel, loci)
    my $runMode = $cgi->param('runMode');

    if($runMode eq 'novel'){
        print STDERR "Novel mode\n";
        print STDERR Dumper($cgi);
    }
    elsif($runMode eq 'pan'){

    }
    elsif($runMode eq 'loci'){

    }
    else{
        print STDERR "Panseq unknown runmode\n";
        exit(1);
    }


}