#!/usr/local/bin/perl

#Cleans Fasta sequences to make sure no errors occur when using Bio:Fasta module
#formats line size and replaces some characters with '_'

package FastaCleaner;

use warnings;
use strict;
use diagnostics;
use Carp;
use IO::File;
use Scalar::Util qw(looks_like_number);

my $fileName;
my $outputName;
my $maxLineSize;

#constructor usage: new FastaCleaner($fileName,$outputName,$maxLineSize);
sub new {
	my $class = shift;
	my $self  = {};

	$fileName    = shift;
	$outputName  = shift;
	$maxLineSize = shift;
	bless $self, $class;
	return $self;
}

#setters
sub setMaxLineSize {
	my $self = shift;
	$maxLineSize = shift;
	unless ( $maxLineSize
		&& ( $maxLineSize == int($maxLineSize) )
		&& ( $maxLineSize > 0 ) )
	{
		confess 'Max Line Size must be an integer greater than 0' . "\n";
	}
}

#getters
sub getMaxLineSize {
	return $maxLineSize;
}

#primary method
sub cleanup {
	my $inputHandle  = new IO::File( '<' . $fileName );
	my $currentLine  = $inputHandle->getline();
	my $outputString = '';
	while ($currentLine) {
		$currentLine =~ s/[\n\r\f]//g;
		if ( $currentLine =~ /^>/ ) {
			$currentLine =~ s/[^>\w\d|.]/_/g;
			$outputString .= $currentLine . "\n";
		}
		else {
			if ( length($currentLine) > $maxLineSize ) {
				my $count = 0;
				my $currentString =
				  substr( $currentLine, $count * $maxLineSize, $maxLineSize );
				while ($currentString) {
					$outputString .= $currentString . "\n";
					$count++;
					if (
						length($currentLine) >
						( ( $count + 1 ) * $maxLineSize ) )
					{
						$currentString =
						  substr( $currentLine, $count * $maxLineSize,
							$maxLineSize );
					}
					else {
						$outputString .= substr(
							$currentLine,
							( $count * $maxLineSize ),
							(
								( $count + 1 ) * $maxLineSize -
								  ( $count * $maxLineSize )
							)
						) . "\n";
						$currentString = '';
					}
				}
			}
			else {
				$outputString .= $currentLine . "\n";
			}
		}
		$currentLine = $inputHandle->getline();
	}
	$inputHandle->close();
	my $outputHandle = new IO::File( '>' . $outputName );
	$outputHandle->print($outputString);
	$outputHandle->close();
}
1;
