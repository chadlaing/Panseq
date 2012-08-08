#!/usr/bin/perl -w

# This module takes in a delta file and creates a simple visualization of the overlapping regions in svg format.
# It creates an embedded svg object within an xml file.
# It is optimized for speed and readability and is reasonably fast.
#
# @author M. Benediktson
# @version 0.1, June 11, 2012
# modified by Chad Laing, June 19, 2012 (object style)

#TODO: Put a scale in the picture. (0 -> 22304304 or sometihng like that)

package Visualization;

use Carp;
use strict;
use warnings;
use FindBin::libs;
use SVG;
use NovelRegion::NovelRegionFinder;
use IO::File;
use Mummer::DeltaBlockFactory;
use FileInteraction::Fasta::SequenceName;
use strict;

sub new {
	my $class = shift;	# Get the class name from the argument array
	my $self = {};
	bless $self,$class;
	$self->_initialize(@_);
	return $self;
}

sub _initialize{
	my $self=shift;
	
	#default values;
	$self->$self->BAR_HEIGHT(25); # This is the height of the bars, in pixels, drawn to the screen
}

sub BAR_HEIGHT{
	my $self=shift;
	$self->{'_BAR_HEIGHT'} = shift // return $self->{'_BAR_HEIGHT'};	
}


# Run the actual program and return an svg file in the form of a string
# Input:
# $filename = the name of the delta file to be read
# (OPTIONAL) $width = the width of the SVG picture
sub run {
	my $self = shift;	# Get the object name from the argument array
	my $filename = shift;	# Get the file name from the argument array

	my $width = do { my $arg = shift; defined($arg) ? $arg : 1000 };	# Width argument is optional

	# Obtain input file information.
	my $FH = IO::File->new('<' . $filename) or die "$!";
	undef $filename;
	my $dbFactory=DeltaBlockFactory->new($FH);

	# Store information in a comparison hash
	my %references = ();
	my %queries = ();
	my $nrf = NovelRegionFinder->new;
	while ( my $block = $dbFactory->nextDeltaBlock(1) ) {    #turns on absolute start/end positions with the true value
		$nrf->updateComparisonHash($block, 'reference' );

		# Add the names to arrays. Perl will not allow duplicate keys, so we don't have to worry about duplicate reference/query names.
		$references{$block->refName} = $block->refLength;	# Store the reference length in the hash for later use
		$queries{$block->queryName} = 1;
	}
	close $FH;	# Close the file. We don't need it anymore.
	undef $dbFactory;
	undef $FH;

	# Construct hash tables of all the references/queries
	# !!SEE THE SUBROUTINE FOR MORE INFO!!
	my %REF = %{&_organize(keys %references)};
	my %QUE = %{&_organize(keys %queries)};
	undef %queries;

	# Find the $max_size value.
	my $max_size = 0;	# The maximum length of the reference.
	for my $ref_key (keys %REF) {
		my $size = 0;
		foreach my $ref_sub (@{$REF{$ref_key}}) {
			$size += $references{$ref_sub};
		}
		if ($size > $max_size) {
			$max_size = $size;
		}
	}

	# Create an SVG object with size $width X $height
	# This is essentially the "canvas" we are going to be drawing to
	my $svg_height = ($self->BAR_HEIGHT * 2 * (keys %QUE)) + $self->BAR_HEIGHT + 1;
	my $svg= SVG->new(width=>$width,height=>$svg_height);
	undef $svg_height;

	# Set the style for all bar objects. It is probably preferable to avoid using the 'stroke' setting
	my %bar = ('fill'=>'black',
		'fill-opacity'=>1);
	my $bar = $svg->group(id=>'bars',style=>\%bar);

	# Set the style for all text objects.
	my %text = ('fill'=>'black',
		'fill-opacity'=>1,
		'text-anchor'=>'middle');
	my $text = $svg->group(id=>'text',style=>\%text);
	my $t1 = $text->text(x=>$width / 2,y=>($self->BAR_HEIGHT / 5) *4); # Put in the center of the screen==

	# Put a Title in the top centre of the image
	# TODO: Do we even need a title for this picture?
	$t1->cdata('Scale: 0 => ' . $max_size);

	# Print titles for the query bars
	my $text_height = 0;
	for my $que_key (keys %QUE) {
		$text_height += $self->BAR_HEIGHT * 2;
		# Print a title for the bar
		$svg->text(x=>$width / 20, y=>$text_height - ($self->BAR_HEIGHT/5), -cdata=>$que_key); #=
	}
	undef $text_height;

	# Do the calculations and draw the queries to the screen. 
	for my $ref_key (keys %REF) {		# For each key in the REF hash.
		my $offset = $width / 20;		# The offset from the left of the screen.
		foreach my $ref_sub (@{$REF{$ref_key}}) {	# for each element of the array at that value of the REF hash.

			my $bar_height = 0;	# used for the bar height.

			for my $que_key ( keys %QUE ) {	# For each key in the QUE hash.

				my @arr = ();	# Initialize an array to hold all the info for the current primary query.
				my $count = 0;	# Initialize a counter for the array.

				# Add all the ranges within a sub-query to a 2D array.
				foreach my $que_sub (@{ $QUE{$que_key} }) {
					# Only continue if this function returns a non-zero value.
					if (my $comp = $nrf->comparisonHash->{$ref_sub}->{$que_sub}) {
						my @ranges = grep (length, split(/\,/, $comp));	# Split and remove all null elements.

						foreach my $token (@ranges) {
							my @temp = split(/\.+/, $token);	# Split on periods.

							# Resize the elements so that they fit in our image.
							$temp[0] = &_resize($temp[0],$max_size,$width);
							$temp[1] = &_resize($temp[1],$max_size,$width);

							push @{$arr[$count]},@temp;	# Add the elements to our array.
							$count++;
						}
					}
				}	
				undef $count;		

				@arr = @{&_merge(@arr)};	# Sort and merge the ranges.

				# Print a bar to the screen.
				$bar_height += ($self->BAR_HEIGHT * 2);
				for (my $i = 0; $arr[$i]; $i++) {
					my $w = $arr[$i][1] - $arr[$i][0];	# Width = End - Start
					if ($w) {				# Ignore zero width sections.
						$bar->rectangle(x=>$offset + $arr[$i][0],	# Left Offset + Start
								y=>$bar_height,
								width=>$w,
								height=>$self->BAR_HEIGHT);
					}
				}
			}
			$offset += &_resize($references{$ref_sub},$max_size,$width);	# Increase the offset from the left of the screen by the length of the current reference.
		}
	}
	my $out = $svg->xmlify;		# Convert to an XML file.
	return $out;			# Preturn a very long string with the xml
}

##=====     BEGIN SUBROUTINE DEFINITIONS     =====##

# Construct a hash of arrays.
# This is a little complicated, so I'll explain it in depth. What we're doing is gathering up all the unique sequences. We are then 
# making each of those unique sequences a key in our hash table. This key points to an array of all sub-sequences (or pieces of 
# sequences) that are part of that unique sequence.
# Input:
# @sequences = an array of all the sub-sequence names.
# Returns:
# A hash reference
sub _organize {
	my (@sequences) = @_;	#Put the arguments into variables.
	my %REAL_ref = ();
	foreach my $value (@sequences) {
		my $seqName = SequenceName->new($value)->name;
		# If the sequence name is not a key in our hash table, Perl will add it and push the value.
		# If the sequence name is already in our hash table, Perl will simply push the value
		push @{ $REAL_ref{$seqName} },$value;
	}
	return \%REAL_ref;
}

# Merge a series of ranges. This will make drawing to the screen much easier.
# This function has O(n) complexity and is about as efficient as we can make it.
# Input: 
# @array = an array of ranges. Ranges are represented as arrays of length 2. 	
# Position 1 = start of range. Position 2 = end of range.
# Returns:
# An array reference
sub _merge {
	my (@array) = @_;	#Put the arguments into variables
	@array = sort {$a->[0] <=> $b->[0]} @array;	# Sort the ranges by first element. 	
							# Example: 	1-5, 8-20, 3-6 would become
							# 		1-5, 3-6, 8-20
	for (my $x = 0; $array[$x+1];) {
		# If the next range overlaps with the first range, merge the two ranges.
		# For example, if we had: x = 1-9 and y = 7-12, we could merge them to create a unified range of 1-12.
		# In addition, if we had: x = 1-9 and y = 9-12, we could merge those to create 1-12.
		# And if we had: x = 1-9 and y = 10-12, we could merge those to create 1-12.
		# Since we are dealing with pixel values, we can safely merge all these without creating errors or discrepancies.
		if ($array[$x][1] + 1 >= $array[$x+1][0]) {
			$array[$x][1] = $array[$x+1][1];	# Merge the two items into the first and then...
			splice(@array, $x+1, 1);		# Delete the second item from the array
		}
		# If the range doesn't overlap, we are done with this set. Increment the counter and move on to the next one.
		else {
			$x++;
		}
	}
	return \@array;
}

# Adjust the size of a value so that it fits in our image properly.
# Input:
# $value = value to be resized
# $max = maximum size of bar in our picture
# $width = the width of the page
# Returns:
# An integer
sub _resize {
	my ($value,$max,$width) = @_;
	# Divide by the maximum bar size.
	$value /= $max;
	# Leave a bit of whitespace on either side of the picture.
	$value *= $width * 0.9;
	# Round to the nearest integer.
	$value = int($value + 0.5);

	return $value;
}
1;
