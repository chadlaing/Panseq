#!/usr/bin/env perl
package Interface::Scripts::Panseq;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use Interface::Scripts::ServerSettings;
use File::Path qw/make_path/;
use IO::File;
use Interface::Scripts::SendMailWrapper;

use parent 'Interface::Scripts::Panseq_Super';

sub setup{
	my $self=shift;


	$self->start_mode('home');
	$self->run_modes(
		'analyses'=>'analyses',
		'submit'=>'submit',
		'home'=>'home',
		'contact'=>'contact'
	);

	$self->serverSettings($self->_loadServerSettings("$FindBin::Bin/../../serverSettings.txt"));
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub serverSettings{
	my $self=shift;
	$self->{'_serverSettings'}=shift // return $self->{'_serverSettings'};
}

sub batchFile{
	my $self=shift;
	$self->{'_batchFile'}=shift // return $self->{'_batchFile'};
}

sub home{
	my $self=shift;

	my $template=$self->load_tmpl('home.tmpl',die_on_bad_params=>0);
	return $template->output();
}

sub analyses{
	my $self=shift;

	my $template = $self->load_tmpl('analyses.tmpl', die_on_bad_params=>0);
	return $template->output();
}

sub contact{
	my $self=shift;

	my $template = $self->load_tmpl('contact.tmpl', die_on_bad_params=>0);
	return $template->output();
}

sub submit{
	my $self=shift;

	my $pid = fork();
    if(!defined $pid){
            die "cannot fork process!\n $!";
    };


    if($pid){
        #display alldone!
        my $LOOP_VAL_REF = $self->_getFormSettings();
		my $template = $self->load_tmpl('submit.tmpl', die_on_bad_params=>0);
		$template->param(LOOP_VAL => $LOOP_VAL_REF);
		return $template->output();
    }
    else{
    	#launch the panseq program
    	close STDERR;
    	close STDOUT;
    	$self->_launchPanseq();
	}
}


sub _launchPanseq{
	my $self=shift;

	$self->serverSettings->baseDirectory($self->_createBaseDirectoryName());
	$self->_createQueryReferenceDirectory();
	my ($q,$batchFile) = $self->_createBatchFile();
	$self->_executePanseqSystemCall($batchFile);
	$self->_sendEmail($q);
}


sub _sendEmail{
	my $self=shift;
	my $q=shift;

	my $mail = Interface::Scripts::SendMailWrapper->new($self->serverSettings->emailDefaultsFile);
	$mail->to($q->param('email'));
	$mail->downloadLink($self->serverSettings->baseDirectory . 'panseq_result.zip');
	$mail->sendTheEmail();
}




sub _executePanseqSystemCall{
	my $self = shift;
	my $batchFile = shift;

	my $systemLine = "perl $FindBin::Bin/../../panseq.pl $batchFile";
	system($systemLine);
	return 1;
}


sub _createQueryReferenceDirectory{
	my $self=shift;

	#all query/reference directories are relative to the fastaFile directory
	#fastaFileDirectory/BaseName/Query
	#fastaFileDirectory/BaseName/Reference

	my $fastaBase = $self->serverSettings->fastaFileDirectory . $self->serverSettings->baseDirectory . '/';
	my $queryDir = $fastaBase . 'query/';
	my $referenceDir = $fastaBase . 'reference/';

	#with File::Path
	#we don't want any permission issues, so don't downgrade the directory permissions
	umask(0); 
	make_path($fastaBase) or die "Couldn't create fastaBase $fastaBase\n";
	make_path($queryDir) or die "Couldn't create query directory at $fastaBase\n";
	make_path($referenceDir) or die "Couldn't create reference directory at $fastaBase\n";

	#set serverSettings variables
	$self->serverSettings->queryDirectory($queryDir);
	$self->serverSettings->referenceDirectory($referenceDir);
}

sub _createStrainSymLink{
	my $self = shift;
	my $type=shift;
	my $strainName = shift;

	#get query/ref directory
	my $linkedLocation; 
	if($type =~ m/^query/){
		$linkedLocation = $self->serverSettings->queryDirectory;
	}
	else{
		$linkedLocation = $self->serverSettings->referenceDirectory;
	}

	my $actualLocation = $self->serverSettings->fastaFileDirectory . $strainName;
	$linkedLocation .= $strainName;

	#create the link
	symlink($actualLocation, $linkedLocation);
}
	

sub _createBaseDirectoryName{
    my $self = shift;
  
    #use random number as well as localtime to ensure no directory overlap
    my $randomNumber = int( rand(8999) ) + 1000;
    my $directory    = localtime . $randomNumber;
    $directory =~ s/[\W]//g; 
  
    return $directory;
}

sub _createBatchFile{
	my $self=shift;

	my $q = $self->query();
	my @params = $q->param();

	my $batchFile = $self->serverSettings->outputDirectory . $self->serverSettings->baseDirectory . '.batch';
	my $batchFH = IO::File->new('>' . $batchFile) or die "Could not create batch file$!\n";

	foreach my $setting(@params){	
		#these are the strains included in the database
		if($setting =~ m/^(querySelected|referenceSelected)/){
		
			#query/ref directory should already exist in serverSettings
			my $stringOfNames = $q->param("$setting");
	
			while($stringOfNames =~ /([\w\.\-]+)/gc){
				my $fileName = $1;
				$self->_createStrainSymLink($setting,$fileName);
			}
			next;
		}


		#these are the strains uploaded by the user
		if($setting =~ m/^(queryFiles|referenceFiles)/){
			$self->_uploadFiles($setting,$q);
		}

		$batchFH->print($setting . "\t" . $q->param($setting) . "\n");
	}

	$batchFH->print('blastDirectory' . "\t" . $self->serverSettings->blastDirectory . "\n");
	$batchFH->print('mummerDirectory' . "\t" . $self->serverSettings->mummerDirectory . "\n");
	$batchFH->print('muscleExecutable' . "\t" . $self->serverSettings->muscleExecutable . "\n");
	$batchFH->print('numberOfCores' . "\t" . $self->serverSettings->numberOfCores . "\n");
	$batchFH->print('baseDirectory' . "\t" . $self->serverSettings->outputDirectory . $self->serverSettings->baseDirectory . "/\n");
	$batchFH->print('queryDirectory' . "\t" . $self->serverSettings->queryDirectory . "\n");
	$batchFH->print('referenceDirectory' . "\t" . $self->serverSettings->referenceDirectory . "\n");
	$batchFH->print('novelRegionFinderMode' . "\tno_duplicates\n"); #hard coded, until interface options exist
	$batchFH->close();

	return ($q,$batchFile);
}

sub _uploadFiles{
	my $self =shift;
	my $type = shift;
	my $q =shift;

	my $directory;
	if($type =~ m/^query/){
            $directory= $self->serverSettings->queryDirectory;
    }
    elsif($type =~ m/^reference/){
            $directory = $self->serverSettings->referenceDirectory;
    }
    else{
            print STDERR "Incorrect type sent to uploadFile\n";
            exit(1);
    }
   
   #upload returns a list of filehandles
   #this requires enctype="multipart/form-data"
    my @FHS = $q->upload($type);

    foreach my $FH(@FHS){
    	$self->_uploadFile($directory,$FH)
    }
}

sub _uploadFile{
        my $self = shift; 
        my $directory=shift;
        my $FH = shift;
      
        my $cleanedFile = $FH;
        $cleanedFile =~ s/\W/_/g;

        if ($FH){
        	#need to get the actual handle, otherwise will try to work on the name, and will fail
        	my $inFH = $FH->handle;

            my $outputHandle = IO::File->new( '>' . $directory . $cleanedFile ) or die "Cannot create $directory$cleanedFile\n";

            #upload file using 1024 byte buffer
            my $buffer;
            my $bytesread = $inFH->read( $buffer, 1024 );
            while ($bytesread) {
                    $outputHandle->print($buffer);
                    $bytesread = $inFH->read( $buffer, 1024 );
            }
            $outputHandle->close();
            return 1;
        }
        else {
        		print STDERR "no upload obj\n";
                return 0;
        }
}

sub _getFormSettings{
	my $self =shift;

	my $q = $self->query();
	my @vals = $q->param();
	my $runMode = $q->param('runMode');

	my @loopVals;
	foreach my $val(@vals){
		unless ($self->_keepCurrentSetting($val, $runMode)){
			next;
		}

		my %settings=(
			'VAL'=>$val // 0,
			'INFO'=>$q->param($val) // 0
		);
		push @loopVals, \%settings;
	}
	my %runMode=(
		'VAL'=>'runMode',
		'INFO'=>$runMode
	);
	push @loopVals, \%runMode;
	return \@loopVals;
}

sub _keepCurrentSetting{
	my $self=shift;
	my $key=shift;
	my $mode=shift;

	my %modeHash=(
		'email'=>1,
		'runMode'=>1
	);

	if($mode eq 'pan' || $mode eq 'novel'){
		$modeHash{'querySelected'}=1;
		$modeHash{'nucB'}=1;
		$modeHash{'nucC'}=1;
		$modeHash{'nucD'}=1;
		$modeHash{'nucG'}=1;
		$modeHash{'nucL'}=1;

		if($mode eq 'pan'){
			$modeHash{'coreGenomeThreshold'}=1;
			$modeHash{'percentIdentityCutoff'}=1;
			$modeHash{'fragmentationSize'}=1;
		}
		else{
			$modeHash{'referenceSelected'}=1;
		}
	}
	elsif($mode eq 'loci'){
		$modeHash{'numberOfLoci'}=1;
		$modeHash{'lociFile'}=1;
	}
	else{
		print STDERR "Unknown Panseq mode: $mode";
		exit(1);
	}


	if($modeHash{$key}){
		return 1;
	}
	else{
		return 0;
	}
}

sub _loadServerSettings{
	my $self=shift;
	my $file = shift;

	my $ss = Interface::Scripts::ServerSettings->new($file);
	return $ss;
}

1;




