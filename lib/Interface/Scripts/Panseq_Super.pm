#!/usr/bin/env perl
package Interface::Scripts::Panseq_Super;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";

use parent 'CGI::Application';
use CGI::Application::Plugin::Config::Simple;
use CGI::Application::Plugin::Redirect;
use CGI::Application::Plugin::Session;

sub cgiapp_init{
	my $self = shift;

	#set paths
	#the template path is set using CGI::Application::Dispatch

	#load config file and get settings


	#Session information
	$self->session_config(DEFAULT_EXPIRY => '+8h');
	
}

1;

