#!/usr/bin/env perl


package Interface::Scripts::Common;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";


#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}


sub _initialize{
	my $self=shift;
}



sub _getPageTop{
	my $self=shift;
my $html = q(

<!DOCTYPE HTML>

<html>
<head>
	<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
	<title>Panseq</title>
	<link href="/panseq.css" rel="stylesheet" type="text/css">
	<script src="jquery-1.9.1.js"></script>
	<script src="panseq.js"></script>
</head>	

);
	return $html;
}


sub _getPageBottom{
	my $self = shift;
my $html = q(

</body>
<footer>
&copy; 2012 Public Health Agency of Canada. All Rights Reserved.
</footer>

</html>
);
	return $html;
}

1;