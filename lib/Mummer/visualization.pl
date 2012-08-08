#!/usr/bin/perl -w

use strict;
use FindBin::libs;
use Mummer::Visualization;

my $svg = Visualization->new();
print $svg->run('/home/michael/multi_ref.delta');
