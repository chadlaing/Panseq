#!/usr/local/bin/perl

#Script written to help deal with and properly organize temporary files
#Unless files are specifed to a path they are redirected to a temp directory
#Reverting back to base directory is performed when end() is called
#JC 4/11/11
package TempFolder;
use warnings;
use strict;
use diagnostics;
use Carp;
use IO::File;
use Cwd 'getcwd';
use File::Path;
use Object::Tiny::RW qw{
  rootFolder
  tempFolder
};

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;

	#generate temp directory name
	$self->rootFolder( getcwd() );
	my $tempFolderName = $self->rootFolder . '/temp';
	$self->tempFolder($tempFolderName);
	my $folderCounter = 0;
	while ( -e $self->tempFolder ) {
		$self->tempFolder( $tempFolderName . $folderCounter );
		$folderCounter++;
	}

	#make temp directories
	mkdir( $self->tempFolder, 0777 ) or mkdir( $tempFolderName . $folderCounter, 0777 );
	chdir $self->tempFolder;
	return $self;
}

#resets everything to as if folder never existed at all
sub end {
	my $self = shift;
	chdir $self->rootFolder;
	$self->_removeFolder();
}

#helper methods
#deletes the folder and everything within it
sub _removeFolder {
	my $self = shift;
	rmtree( $self->tempFolder );
}

#Module Testing purposes only
#changes dir but retains files
sub endNoDelete {
	my $self = shift;
	chdir $self->rootFolder;
}
1;
