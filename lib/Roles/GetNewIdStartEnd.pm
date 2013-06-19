#!/usr/bin/env perl

package Roles::GetNewIdStartEnd;

use strict;
use warnings;
use FindBin;
use IO::File;
use lib "$FindBin::Bin/../";
use Role::Tiny;

=head2 _getNewIdStartEnd

Checks for a previous substring ending.
If this exists, updates the new coords relative to the initial sequence,
replacing the old substring ending with the new one.
If not, return the current substring coords and id.

=cut

sub _getNewIdStartEnd{
	my $self = shift;
	my $currId = shift;
	my $currStart =shift;
	my $currEnd = shift;

	#check for previous substring
	if($currId =~ m/(_\((\d+)\.\.(\d+)\))$/){
		my $oldLocation = $1;
		my $oldStart = $2;
		my $oldEnd= $3;

		my $delta = $currEnd - $currStart;
		my $newStart = $oldStart + $currStart -1;
		my $newEnd = $newStart + $delta;

		my $newId = $currId;

		$newId =~ s/\Q$oldLocation\E//;

		return($newId,$newStart,$newEnd);
	}
	else{
		return($currId, $currStart, $currEnd);
	}
}

1;
