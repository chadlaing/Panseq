#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Bio::SeqIO;
use Modules::Fasta::SequenceName;

=pod

=head1 NAME

create_strain_template.pl - Given a directory name, create an HTML::Template file of all strains.

=head1 SYNOPSIS

	perl create_strain_template.pl /my/directory/ > my_output_file.tmpl

=head1 DESCRIPTION

This script generates a static template file listing all strains.
Can be set up to run as a CRON job to keep list current.

=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: 

=head1 AUTHOR

Your name (chadlaing@gmail.com)

=cut

my $directory = $ARGV[0];

opendir( DIRECTORY, $directory ) or die "cannot open directory $directory $!\n";
my @dir = readdir DIRECTORY;
closedir DIRECTORY;

foreach my $file(@dir){
	next if substr( $file, 0, 1 ) eq '.';

	my $fullPathFile = $directory.$file;
	my $inFH =Bio::SeqIO->new(-file=>$fullPathFile, -format=>'fasta') or die "Could not open file $fullPathFile\n";

	my $seq = $inFH->next_seq();
	my $fullHeader = $seq->id() . ' ' . $seq->desc();
	print _createTemplateLine($fullHeader,$file);
	$inFH->close();
}

sub _createTemplateLine{
	my $header =shift;
	my $file = shift;

	#remove any contig reference from the header
	if($header =~ m/\s(\S*(cont|ctg)\S*)\s/i){
		my $contig = $1;
		$header =~ s/$contig//;
	}	

	my $humanLine=$header;
	if($header =~ m/^\S+\s+(.+)/){
		$humanLine =$1;
	}

	return(qq{<li name="$file">$humanLine</li>\n});
}