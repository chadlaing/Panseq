#!/usr/bin/perl
use IO::File;

use strict;
use warnings;

my $treeFile = $ARGV[0];
my $infoFile = $ARGV[1];

unless(defined $infoFile && defined $treeFile){
	print "Missing parameters. Correct usage is:\n" ,
		'perl treeNumberToName.pl <tree_file_name> <info_file_name> > <output_file_name>' . "\n";
	exit(1);
}

my $treeFH = IO::File->new('<' . $treeFile) or die "$treeFile does not exist.\n$!";
my $infoFH = IO::File->new('<' . $infoFile) or die "$!";

#instead of File::Slurp
my $treeAsLine;
while(my $line = $treeFH->getline){
	$treeAsLine .= $line;
}

my $counter=0;
while(my $line = $infoFH->getline){
	if($counter==0){
		$counter++;
		next;
	}
	
	$line =~s/\R//g;
	my @la= split('\t',$line);
	my $number = $la[0];
	my $name = $la[1];	
	
	#error proof name
	$name =~ s/[\s\:]/_/g;
	#print STDERR "number: $number\nname: $name\n";

	#substitute
	if($treeAsLine =~ m/[\)\(\,\s]\Q$number\E[\:\)]/){
		print STDERR "matched\n";
		$treeAsLine =~ s/([\)\(\,\s])\Q$number\E([\:\)])/$1$name$2/;
	}
}

print $treeAsLine;
$treeFH->close();
$infoFH->close();