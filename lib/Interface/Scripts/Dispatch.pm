#!/usr/bin/env perl
package Interface::Scripts::Dispatch;

# mod_rewrite alters the PATH_INFO by turning it into a file system path,
# so we repair it.
#from https://metacpan.org/module/CGI::Application::Dispatch#DISPATCH-TABLE

# $ENV{PATH_INFO} =~ s/^$ENV{DOCUMENT_ROOT}// if defined $ENV{PATH_INFO};

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use parent qw/CGI::Application::Dispatch/;

sub dispatch_args {
    return {
        prefix  => 'Interface::Scripts',
        args_to_new=>{
        	TMPL_PATH=>"$FindBin::Bin/../Templates/"
        },
        table   => [
        	''	=>	{app=>'Panseq',rm=>'home'},
        	'/'	=>	{app=>'Panseq',rm=>'home'},
        	'/home/' => {app=>'Panseq',rm=>'home'},
            '/home' => {app=>'Panseq',rm=>'home'},
            '/analyses/' => {app=>'Panseq',rm=>'analyses'},
            '/analyses' => {app=>'Panseq',rm=>'analyses'},
            '/contact/' => {app=>'Panseq',rm=>'contact'},
            '/contact' => {app=>'Panseq',rm=>'contact'},
            '/submit/'=> {app=>'Panseq', rm=>'submit'},
            '/waiting' => {app=>'Panseq', rm=>'waiting'},
            '/download/:analysis_id' => {app=>'Panseq', rm=>'download'},
            'submit'=>{app=>'Panseq',rm=>'submit'},
            ':app/:rm'        => { },
            'admin/:app/:rm'  => { prefix   => 'MyApp::Admin' },
        ],
    };
}

1;

