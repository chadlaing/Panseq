#!/usr/bin/perl

=pod

=head1 NAME

Modules::Setup::PanseqFiles - Stores all filenames from query/reference directories.
Also provides functionality for combining all query/reference files into a single file.

=head1 SYNOPSIS

	use FindBin::libs;
	use Modules::Setup::PanseqFiles;
	
	my $files = Modules::Setup::PanseqFiles->new(
		'queryDirectory'=>$settings->queryDirectory,
		'referenceDirectory'=>$settings->referenceDirectory
	);

=head1 DESCRIPTION

This module allows for the management of files in the Panseq query and reference folders.
An undefined value indicates there are no files in the folder.

=cut

=head1 ACKNOWLEDGEMENTS

Thanks.

=head1 COPYRIGHT

This work is released under the GNU General Public License v3  http://www.gnu.org/licenses/gpl.html

=head1 AVAILABILITY

The most recent version of the code may be found at: https://github.com/chadlaing/Panseq

=head1 AUTHOR

Chad Laing (chadlaing gmail com)

=head2 Methods

=cut

package Modules::Setup::PanseqFiles;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use File::Copy;
use Log::Log4perl;
use Role::Tiny::With;

with 'Roles::CombineFilesIntoSingleFile';

#object creation
sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}


=head3 queryFileNames

An array ref of all the query files in the query directory.

=cut

sub queryFileNames{
	my $self=shift;
	$self->{'_queryFiles'} = shift // return $self->{'_queryFiles'};
}


=head3 referenceFileNames

An array ref of all the reference files in the reference directory

=cut

sub referenceFileNames{
	my $self=shift;
	$self->{'_referenceFiles'} = shift // return $self->{'_referenceFiles'};
}


=head3 _initialize

Initializes the logger.
Takes in any of the two parameters 'queryDirectory' and 'referenceDirectory'.
Any other option will result in a program halt.
If either options are present, the _gatherFileNames(directory) sub is called (via Roles::CombineFilesIntoSingleFile)
which retrieves the file names from the appropriate directory and stores them as queryDirectory
or referenceDirectory, and the names of the files in each as queryFileNames or referenceFileNames

=cut

sub _initialize{
	my($self)=shift;

    #logging
    $self->logger(Log::Log4perl->get_logger()); 

    $self->logger->debug("Logger initialized in Modules::Setup::PanseqFiles");  

    my %params = @_;

    foreach my $key(keys %params){
    	if($key eq 'queryDirectory'){
    		$self->queryFileNames($self->_getFileNamesFromDirectory($params{$key}));
    		$self->logger->debug("Gathering query file names from $params{$key}");
    	}
    	elsif($key eq 'referenceDirectory'){
    		$self->referenceFileNames($self->_getFileNamesFromDirectory($params{$key}));
    		$self->logger->debug("Gathering reference file names from $params{$key}");
    	}
    	else{
    		$self->logger->info("Param $key is not valid in Modules::Setup::PanseqFiles");
    		exit(1);
    	}
    }


    #the reference files can be left blank if building a pan-genome
    #if this is left undefined, an error will be thrown
    #set to an anon array to allow blank reference, but not blank query files
    unless(defined $self->referenceFileNames){
    	$self->referenceFileNames([]);
    }

}

=head3 singleQueryFile

All files from the query directory are combined into a single fasta file.
If this file exists, the name is returned.
If it does not, the files are combined and the name of the combined file is returned.
Use of Roles::CombineFilesIntoSingleFile

=cut

sub singleQueryFile{
	my $self=shift;

	if(defined $self->_singleQueryFileName){
		return $self->_singleQueryFileName;
	}
	else{
		my $combinedFileName = shift // undef;
		$self->_singleQueryFileName($self->_combineAndSanitizeFastaFiles($self->queryFileNames, $combinedFileName));
	}

}

=head3 singleReferenceFile

All files from the reference directory are combined into a single fasta file.
If this file exists, the name is returned.
If it does not, the files are combined and the name of the combined file is returned.
Use of Roles::CombineFilesIntoSingleFile

=cut

sub singleReferenceFile{
	my $self=shift;

    if(defined $self->_singleReferenceFileName){
        return $self->_singleReferenceFileName;
    }
    else{
        my $combinedFileName = shift // undef;
        $self->_singleReferenceFileName($self->_combineAndSanitizeFastaFiles($self->referenceFileNames, $combinedFileName));
    }

}


=head3 allQueryAndReferenceFilesAsSingleFile

The combination of the query and reference directory into a single file.
If this file exists, the name is returned.
Uses Roles::CombineFilesIntoSingleFile.

=cut

sub allQueryAndReferenceFilesAsSingleFile{
    my $self=shift;

    if(defined $self->_allQueryAndReferenceFilesAsSingleFile){
        return $self->_allQueryAndReferenceFilesAsSingleFile;
    }
    else{
        my $combinedFileName = shift // undef;
        my $filesToCombine = [@{$self->queryFileNames},@{$self->referenceFileNames}];
        $self->_singleReferenceFileName($self->_combineAndSanitizeFastaFiles($filesToCombine, $combinedFileName));
    }
}


sub _allQueryAndReferenceFilesAsSingleFile{
    my $self=shift;
    $self->{'__allQueryAndReferenceFilesAsSingleFile'}=shift // return $self->{'__allQueryAndReferenceFilesAsSingleFile'};
}


=head3 _singleReferenceFileName

Stores the filename of the combined reference files.
Accessed by singleReferenceFile only.

=cut

sub _singleReferenceFileName{
	my $self=shift;
	$self->{'__singleReferenceFileName'} = shift // return $self->{'__singleReferenceFileName'};
}

=head3 _singleReferenceFileName

Stores the filename of the combined query files.
Accessed by singleQueryFile only.

=cut

sub _singleQueryFileName{
	my $self=shift;
	$self->{'__singleQueryFileName'} = shift // return $self->{'__singleQueryFileName'};
}


1;



