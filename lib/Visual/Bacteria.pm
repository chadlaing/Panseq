=head1 NAME

Visual::Bacteria

=head1 SYNOPSIS

	use Visual::Bacteria
	my $object = Visual::Bacteria->new();
	print $object->run('input_filename', 'OPTIONS');

=head1 DESCRIPTION

This module takes in a tab-delimited text file and creates a simple visualization of the factors present
It creates an embedded svg object within an xml file.
The user may select a number of different display options.

=head2 METHODS

=over 12

=item C<new>

Returns a new Visual::Bacteria object.

=item C<run>

Run the actual program and return an svg file in the form of a string
Input:
$filename = the name of the delta file to be read
Optional (and in any order after the mandatory arguments):
$colour = Either 'red', 'green' or 'blue'. Defaults to 'blue'.
$shape = Either 'circle' or 'square'. Defaults to 'circle'.
$size = Either 'small', 'medium' or 'large'. Defaults to 'medium'.
In addition, the order of the header can be set by entering a filename. This must be a ".txt" file and the names must be newline delimited.

=item C<_check_input>

This subroutine is used to check for errors in the input file.
Upon reaching an error, the program will revert to the default
sorting settings.

NOTE: It might be a good idea to display specific error messages for the errors.

INPUT:
$input = reference to the original array of names
$ouput = reference to the sorted array of names
RETURNS:
This subroutine modifies both the $input and the $output arrays

=item C<_grid>

Prints a grid of an appropriate size to the scre$hash, $keys, $order, $size, $colours, $shape, $writeren.

INPUT:
$horizontal = the number of horizontal lines to draw
$vertical = the number of vertical lines to draw
$size = the size of the grid boxes
$writer = the XML::Writer object
RETURNS:
nothing

=item C<_shapes>

This subroutine draws the shapes onto the chart at their correct positions.
It also prints the appropriate mouseover text.

INPUT:
$hash = a reference to a hash of genomes
$keys = a reference to an array of sorted genome keys
$order = a reference to an array with the current ordering of the headers
$width = the width of the page
$size = the size of the shapes
$colours = a reference to a 2d array of colours
$shape = the currently selected shape
$writer = the XML::Writer object 
RETURNS:
nothing

=item C<_header>

This subroutine draws the header text.

INPUT:
$names = a reference to an array of sorted names
$size = the size of the shapes
$writer = the XML::Writer object 
RETURNS:
nothing

=item C<_header_mousover>

This subroutine draws the header text mouseover.
NOTE: We have to do this in a seperate loop AFTER drawing the squares, or the text won't be readable.

INPUT:
$names = a reference to an array of sorted names
$width = the width of the page
$size = the size of the shapes
$writer = the XML::Writer object 
RETURNS:
nothing

=head1 AUTHOR

M. Benediktson

=head1 VERSION

0.2, July 9, 2012

=cut

#!/usr/bin/perl -w

package Visual::Bacteria;

use FindBin;
use lib "$FindBin::Bin/../../lib";
use Carp;
use XML::Writer;
use XML::Writer::String;
use IO::File;
use strict;

use constant BOX_SIZE => 25;		# This is the height of the bars, in pixels, drawn to the screen.
use constant CHART_MARGIN => 165;	# This is the space allowance for the text that would otherwise be drawn offscreen.

use constant COLOUR_LINES => '#777777'; # This is the colour of the grid lines


sub new {
	my $class = shift;	# Get the class name from the argument array
	my $self = {};
	bless $self,$class;
	return $self;
}

# Run the actual program and return an svg file in the form of a string
# Input:
# $filename = the name of the delta file to be read
# Optional (and in any order after the mandatory arguments):
# $colour = Either 'red', 'green' or 'blue'. Defaults to 'blue'.
# $shape = Either 'circle' or 'square'. Defaults to 'circle'.
# $size = Either 'small', 'medium' or 'large'. Defaults to 'medium'.
# In addition, the order of the header can be set by entering a filename. This must be a ".txt" file and the names must be newline delimited.
sub run {
	my $self = shift;	# Get the object name from the argument array
	my $filename = shift;	# Get the file name from the argument array
	my $shape = 'circle';	# Set the default shape for the symbols
	my $size = BOX_SIZE;	# Set the defualt size for the chart
	my @sorted_names = ();	# Set the ordering array

	# These are all the colours we are going to be using in order of: "Primary Colour", "Stroke Colour", "Highlight Colour"
	my @colours = (['#006FFF','#3299CC','#3232CD'],		# Blues
		       ['#A62A2A','#FF2400','#8C1717'],		# Reds
		       ['#32CD32','#00FF7F','#238E23'],		# Greens
		       ['#CC3299','#FF6EC7','#8E236B'],		# Pinks
		       ['#6B4226','#855E42','#5C3317'],		# Browns
		       ['#D9D919','#FFFF00','#EAEAAE'],		# Yellows
		       ['#660099','#9933CC','#330066'],		# Purples
		       ['#FF5721','#FF9955','#D43D1A']);	# Oranges

	# Read in all the optional arguments
	# We can read in the arguments in any order to allow the user some greater flexibility.
	while (defined(my $arg = shift)) {
		# Shape of the symbols
		if ($arg eq 'circle' || $arg eq 'square') {
			$shape = $arg;
		}
		# Size of the entire chart
		elsif ($arg eq 'small' || $arg eq 'medium' || $arg eq 'large') {
			if ($arg eq 'small') {
				$size -= 5;
			}
			elsif ($arg eq 'large') {
				$size += 5;
			}
		}
		# Colour of the first symbol that appears on the chart. (This is the colour that '1' will be coloured.)
		elsif ($arg eq 'blue' || $arg eq 'red' || $arg eq 'green') {
			if ($arg eq 'red') {
				@colours[0,1] = @colours[1,0];
			}
			elsif ($arg eq 'green') {
				@colours[0,2] = @colours[2,0];
			}
		}
		# Open a file with an ordering for the header
		elsif ($arg =~ /\.txt$/i){
			my $inFH = IO::File->new('<'. $arg) or die "$!";			

			while (my $line = $inFH->getline) {
				chomp($line);
				push(@sorted_names,$line);
			}
			$inFH->close();
		}
		else {
			#Invalid or unreadable argument
		}
	}

	# These variables are used throughout most of the program
	my @names = ();
	my %genome = ();

	#open(INPUT_FILE, '<' . $filename) or die "Cannot open $filename";
	my $inFH = IO::File->new('<' . $filename) or die "Cannot open $filename";
	
	# Construct a hash of string values. The string values are essentially binary values and all of equal length. ( 0 or 1 only)
	my $count=0;
	while(my $line = $inFH->getline) {
	 	chomp($line);

		my @arr = grep(length, split(' ', $line));	# Split at "tab" and remove all null elements.

		# Special case for the first line. We need to get the names
		if ($count == 0) {
			@names = @arr;
			$count++;
		}
		else {
			for (my $i = 1; $i < @arr; $i++) {
				$genome{$arr[0]} .= $arr[$i];	# Add the elements to the string.
			}
		}
	}
	undef $filename;
	$inFH->close();

	# This block is used to check for errors in the input text file. (If one is given. If not, just sort the array of @names.)
	# Errors in the text file will cause the program to revert to the default, sorted setting.
	if (@sorted_names > 0) {
		$self->_check_input(\@names, \@sorted_names);
	}
	else {
		@sorted_names = sort(@names);
	}

	# We need to know how the names array is sorted. If we know how it's sorted, we can sort the headers.
	# So this is where we create an 'ordering' array that tells us how the names array is sorted.
	my @ordering;
	my %index;
	@index{@names} = (0..$#names);	
	for (my $i = 0; $i < @names; $i++) {
		push(@ordering, $index{$sorted_names[$i]});
	}
	undef @names;	# We use @sorted_names below so there's no sense in keeping another array with the same items.
	undef %index;

	my @genome_keys = sort(keys %genome);	# Make a list of all the keys and sort them alphabetically.

	# Create an SVG object with size $width X $height
	# This is essentially the "canvas" we are going to be drawing to.
	my $width = ($size * (scalar @sorted_names)) + CHART_MARGIN + 100;	# Width of the entire SVG. Used later on.
	my $s = XML::Writer::String->new();				# The string that we are printing the file to.
	my $writer = XML::Writer->new( OUTPUT => $s);

	$writer->xmlDecl('UTF-8');
	$writer->doctype('svg', '-//W3C//DTD SVG 20001102//EN','http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd');
	$writer->startTag('svg',
		           height => ($size * (@genome_keys + 2)) + CHART_MARGIN,
		           width  => $width,
			   xmlns  => 'http://www.w3.org/2000/svg');

	# Print the grid to the screen.
	$self->_grid(scalar(keys %genome), scalar(@sorted_names), $size, $writer);

	$self->_header(\@sorted_names, $size, $writer);	# Draw the header text
	$self->_shapes(\%genome, \@genome_keys, \@ordering, $width, $size, \@colours, $shape, $writer);	# Draw the shapes on the chart
	$self->_header_mouseover(\@sorted_names, $width, $size, $writer);	# Draw the header text mouseover

	$writer->endTag('svg');
	$writer->end();

	return $s->value();
}

##=====     BEGIN SUBROUTINE DEFINITIONS     =====##


# This subroutine is used to check for errors in the input file.
# Upon reaching an error, the program will revert to the default
# sorting settings.
#
# NOTE: It might be a good idea to display specific error messages for the errors.
#
# INPUT:
# $input = reference to the original array of names
# $ouput = reference to the sorted array of names
# RETURNS:
# This subroutine modifies both the $input and the $output arrays
sub _check_input {
	my $self=shift;
	my ($input, $output) = @_;		# Get the arguments.

	# Check if they're even the same size. If they're not, revert to the default sorting method.
	if (scalar(@$output) != scalar(@$input)) {
		@$output = sort(@$input);
	}

	my $error = 0;
	my @input_copy = @$input;	# Make a copy of the input array. We'll need to use this later.

	# Make sure that all values in the output array exist in the input array
	for (my $i = 0; $i < @$input && $error == 0; $i++) {
		for (my $j = 0, $error = 1; $j < @$output; $j++) {
			# If we find a match...
			if ($output->[$j] eq $input_copy[$i]) {
				$error = 0;			# remove the error flag,
				$input_copy[$i] = undef;	# delete the original (This makes sure we don't have any copies...),
				last;				# and exit the loop.
			}
		}
	}
	# If an error exits, revert to the default sorting;
	if ($error == 1) {
		@$output = sort(@{$input});
	}
	undef $error;
}

# Prints a grid of an appropriate size to the scre$hash, $keys, $order, $size, $colours, $shape, $writeren.
#
# INPUT:
# $horizontal = the number of horizontal lines to draw
# $vertical = the number of vertical lines to draw
# $size = the size of the grid boxes
# $writer = the XML::Writer object
# RETURNS:
# nothing
sub _grid {
	my $self=shift;
	my ($horizontal, $vertical, $size, $writer) = @_;	# Get the arguments.

	# Print horizontal guidelines.
	for (my $i = 0; $i < $horizontal; $i++) {
		$writer->startTag('rect', 
				  'x'=>($size*0.5),
				  'y'=>($size*1.5)+$i*$size+CHART_MARGIN,
				  'fill'=>COLOUR_LINES, 
				  'width'=>$size*($vertical-1),
				  'height'=>1);
		$writer->endTag('rect');		
	}

	# Print vertical guidelines.
	for (my $i = 0; $i < $vertical; $i++) {
		$writer->startTag('rect', 
				  'x'=>($size*0.5)+($i*$size),
				  'y'=>($size*1.5)+CHART_MARGIN,
				  'fill'=>COLOUR_LINES, 
				  'width'=>1,
				  'height'=>$size*($horizontal-1));
		$writer->endTag('rect');		
	}
}

# This subroutine draws the shapes onto the chart at their correct positions.
# It also prints the appropriate mouseover text.
#
# INPUT:
# $hash = a reference to a hash of genomes
# $keys = a reference to an array of sorted genome keys
# $order = a reference to an array with the current ordering of the headers
# $width = the width of the page
# $size = the size of the shapes
# $colours = a reference to a 2d array of colours
# $shape = the currently selected shape
# $writer = the XML::Writer object 
# RETURNS:
# nothing
sub _shapes {
	my $self=shift;
	my ($hash, $keys, $order, $width, $size, $colours, $shape, $writer) = @_;	# Get the arguments

	# Draw the shapes.
	for (my $i = 0, my $height = $size; $i < @$keys; $i++, $height += $size) {
		my $str = $hash->{$keys->[$i]};
		for (my $j = 0; $j < length($str); $j++) {
			# Check if the value is a one. If it is, draw the square.
			my $type = substr($str, $order->[$j], 1);
			if ($type) {

				# Set the shape of the objects
				my $radius = 0;
				if ($shape eq 'circle') {
					$radius = $size;
				}

				my $colour1 = $colours->[0][0];
				my $colour2 = $colours->[0][1];
				my $colour3 = $colours->[0][2];
				# Make sure it's within our colour range	
				if ($type <= @$colours) {
					$colour1 = $colours->[$type-1][0];
					$colour2 = $colours->[$type-1][1];
					$colour3 = $colours->[$type-1][2];
				}			

				# REMEMBER!! The id chosen here must be used later on when we do the mouseover.
				$writer->startTag('rect',
						  'id'=>'item' . $i . ':' . $j, 
						  'x'=>$j * $size + 1,
						  'y'=>CHART_MARGIN + $height + 1,
						  'width'=>$size-1,
						  'height'=>$size-1,
						  'fill'=>$colour1,
						  'rx'=>$radius,
						  'ry'=>$radius,
						  'opacity'=>0.9,
						  'stroke'=>$colour2,
						  'stroke-width'=>2);
				# Change the colour of the box when we mouse over it.
				$writer->startTag('set',
						  'attributeName'=>'fill', 
						  'to'=>$colour3,
						  'begin'=>'mouseover', 
						  'end'=>'mouseout');
				$writer->endTag('set');
				$writer->endTag('rect');
			}
		}
	}

	# Draw the mouseover text for each shape.
	# NOTE: We have to do this in a seperate loop AFTER drawing the shapes, or the text won't be readable.
	for (my $i = 0, my $height = $size; $i < @$keys; $i++, $height += $size) {
		my $str = $hash->{$keys->[$i]};
		for (my $j = 0; $j < length($str); $j++) {
			# Check if the value is a one. If it is, draw the text.
			if (substr($str, $order->[$j], 1)) {
				# Make the invisible text for the shape.
				my $x = $j * $size;
				my $y = $height - ($size/5) + CHART_MARGIN;

				# Make sure the text doesn't go offscreen
				if ($x > $width/2) {
					$x -= 550 + $size;
				}

				# Create an embedded html object for text
				$writer->startTag('foreignObject', 
						  'x'=>$x + $size, 
						  'y'=>$y - 12, 
						  'width'=>550, 
						  'height'=>300, 
						  'visibility'=>'hidden');
				# Draw the html text
				$writer->startTag('p', 
						  'xmlns'=>'http://www.w3.org/1999/xhtml', 
						  'style'=>'background-color:white; 
							    border:1px solid;
							    border-radius:3px;
							    word-wrap:break-word');
				$writer->characters($keys->[$i]);
				$writer->endTag('p');

				# Create the mouseover effect for the text
				$writer->startTag('set',
						  'attributeName'=>'visibility', 
						  'to'=>'visible', 
						  'begin'=>'item'. $i . ':' . $j . '.mouseover', 
						  'end'=>'item' . $i . ':' . $j . '.mouseout');
				$writer->endTag('set');

				$writer->endTag('foreignObject');

			}
		}
	}
}

# This subroutine draws the header text.
#
# INPUT:
# $names = a reference to an array of sorted names
# $size = the size of the shapes
# $writer = the XML::Writer object 
# RETURNS:
# nothing
sub _header {
	my $self=shift;
	my ($names, $size, $writer) = @_;

	# Print the header text.
	for (my $i = 0; $i < @$names; $i++) {
		my $x = $i * $size + 12;	# Add a bit of an offset here to center the text
		my $y = ($size*0.75) + CHART_MARGIN;
		$writer->startTag('text',
				  'id'=>'text' . $i,
				  'font-size'=>16,
				  'transform'=>'rotate(-45 ' . $x . ',' . $y . ')',
				  'x'=>$x,
				  'y'=>$y);
		$writer->characters(substr($names->[$i], 0, 20));	# Print only the first few characters

		# If the string is longer, show that it continues with ellipses
		# Note that we have to compare ASCII values here. Otherwise '0's mess things up.
		# if (ord(substr($names->[$i], 20, 1)) != 0) {
		# 	$writer->characters('...');
		# }

		#Chad Laing edit August 09, 2012
		#test whether there are greater than 20 chars in string
		#if so, add the ... to the name
		if(length($names->[$i]) > 20){
			$writer->characters('...');
		}
		$writer->endTag('text');		
	}
}

# This subroutine draws the header text mouseover.
# NOTE: We have to do this in a seperate loop AFTER drawing the squares, or the text won't be readable.
#
# INPUT:
# $names = a reference to an array of sorted names
# $width = the width of the page
# $size = the size of the shapes
# $writer = the XML::Writer object 
# RETURNS:
# nothing
sub _header_mouseover {
	my $self=shift;
	my ($names, $width, $size, $writer) = @_;	# Get the arguments

	for (my $i = 0; $i < @$names; $i++) {
		# We don't need the mouseover text if it's short.
		#Chad Laing edit August 09, 2012, testing length of string
		if (length($names->[$i]) > 20) {
			my $x = $i * $size + 12;
			my $y = ($size*0.75) + CHART_MARGIN;

			# Make sure the text doesn't go offscreen
			if ($x > $width/2) {
				$x -= 550 + $size;
			}

			# Create an embedded html object for text
			$writer->startTag('foreignObject', 
					  'x'=>$x + $size, 
					  'y'=>$y - 12, 
					  'width'=>550, 
					  'height'=>300, 
					  'visibility'=>'hidden');
			# Draw the html text
			$writer->startTag('p', 
					  'xmlns'=>'http://www.w3.org/1999/xhtml', 
					  'style'=>'background-color:white; 
						    border:1px solid;
						    border-radius:3px;
						    word-wrap:break-word');
			$writer->characters($names->[$i]);
			$writer->endTag('p');

			# Create the mouseover effect for the text
			$writer->startTag('set',
					  'attributeName'=>'visibility', 
					  'to'=>'visible', 
					  'begin'=>'text'. $i . '.mouseover', 
					  'end'=>'text' . $i . '.mouseout');
			$writer->endTag('set');

			$writer->endTag('foreignObject');
		} 		
	}
}
1;
