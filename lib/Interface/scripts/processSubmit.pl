#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use CGI;
use Data::Dumper;
use File::Path qw/make_path/;
use Net::FTP;
use Archive::Extract;

my $cgi = CGI->new();
my $NCBI_PREFIX = '/genomes/all/';

my $pid = fork();
if(!defined $pid){
    die "cannot fork process!\n $!";
};

if($pid){
    #the redirect url goes here
    print $cgi->redirect("/page/index.html");

}
else{


    print STDERR "processing the submission\n";

    #we need to determine what run mode (pan, novel, loci)
    my $runMode = $cgi->param('runMode');
    my $serverSettings = _loadServerSettings();

    if($runMode eq 'novel'){
        foreach my $k(keys %{$cgi->{'param'}}){
            print STDERR "$k\n";
        }


        #for a novel run we need both query / reference directories
        my $newDir = $serverSettings->{'outputDirectory'} . _createBaseDirectoryName();
        my $queryDir = $newDir . 'query/';
        my $refDir = $newDir . 'reference/';
        my $resultsDir = $newDir . 'results/';

        _createDirectory($newDir);
        _createDirectory($queryDir);
        _createDirectory($refDir);

        #create a hash of settings
        #base directory is used in the Panseq program to mean results dir
        #in this hash the results dir is one level down from the base dir
        #baseDirectory in Panseq really should have been called results,
        #but it is too late for that change now

        my %runSettings = (
            queryDirectory => $queryDir,
            referenceDirectory => $refDir,
            mummerDirectory => $serverSettings->{'mummerDirectory'},
            blastDirectory => $serverSettings->{'blastDirectory'},
            numberOfCores => $serverSettings->{'numberOfCores'},
            baseDirectory => $resultsDir,
            numberOfCores => $serverSettings->{'numberOfCores'},
            muscleExecutable => $serverSettings->{'muscleExecutable'},
            outputDirectory => $newDir
        );

        my $batchFile = _createBatchFile(\%runSettings);
        _downloadUserSelections(\%runSettings);
        _runPanseq($serverSettings, $batchFile);

    }
    elsif($runMode eq 'pan'){

    }
    elsif($runMode eq 'loci'){

    }
    else{
        print STDERR "Panseq unknown runmode\n";
        exit(1);
    }
}


sub _runPanseq{
    my $serverSettings = shift;
    my $configFile = shift;
    my $systemLine = 'perl ' . $serverSettings->{'panseqDirectory'} . 'panseq.pl ' . $configFile;
    my $systemExit = system($systemLine);

    if($systemExit == 0){
        #panseq finished properly, send to results html

    }
    else{
        #error in system call
        #send to error html
    }
}



sub _downloadUserSelections{
    my $runSettings = shift;

    my @ncbiQueryGenomes;
    my @ncbiReferenceGenomes;
    foreach my $p(keys %{$cgi->{'param'}}){
        if($p =~ m/^q_(.+)/){
            push @ncbiQueryGenomes, $1;
        }
        elsif($p =~ m/^r_(.+)/){
            push @ncbiReferenceGenomes, $1;
        }
        else {
            print STDERR "$p no match\n";
        }
    }


    #sets up the parameters for ncbi ftp connection
    my $host = 'ftp.ncbi.nlm.nih.gov';
    #constructs the connection
    my $ftp = Net::FTP->new($host, Debug => 1,Passive => 1) or die "Cannot connect to genbank: $@";
    #log in as anonymous, use email as password
    $ftp->login("anonymous",'chadlaing@inoutbox.com') or die "Cannot login ", $ftp->message;


    $ftp->binary();
    my @allDownloadedFiles;
    foreach my $q(@ncbiQueryGenomes){
        my $ncbiFile = $NCBI_PREFIX . $q . '/' . $q . '_genomic.fna.gz';
        my $localFile = $runSettings->{'queryDirectory'} . $q;

        push @allDownloadedFiles, $localFile;
        $ftp->get($ncbiFile, $localFile) or die "Cannot get $q", $ftp->message;
    }

    foreach my $r(@ncbiReferenceGenomes){
        my $ncbiFile = $NCBI_PREFIX . $r . '/' . $r . '_genomic.fna.gz';
        my $localFile = $runSettings->{'referenceDirectory'} . $r;

        push @allDownloadedFiles, $localFile;
        $ftp->get($ncbiFile, $localFile) or die ("Cannot get $r" . $ftp->message);
    }
    $ftp->ascii();


    #extract them all
    foreach my $f(@allDownloadedFiles){
        my $extracter = Archive::Extract->new('archive'=>$f, type=>'gz');
        $extracter->extract(to=>$f . '.fna') or die ("Could not extract $f\n" and next);
        unlink $f;
    }

}


sub _createBatchFile{
    my $paramRef = shift;

    my $batchFile = $paramRef->{'outputDirectory'} . 'panseq.batch';
    open(my $batchFH, '>', $batchFile) or die "Could not create $batchFile\n$!";

    foreach my $k(keys %{$paramRef}){
        $batchFH->print($k . "\t" . $paramRef->{$k} . "\n");
    }

    $batchFH->close();
    return $batchFile;
}



sub _loadServerSettings{
    my $symlinkFile = './server.conf';

    my %settings;
    open(my $inFH, '<', $symlinkFile) or die "$!\n";

    while(my $line = $inFH->getline){
        $line =~ s/\R//g;
        my @la = split(/\s+/, $line);

        $settings{$la[0]}=$la[1];
    }
    return \%settings;
}



sub _createBaseDirectoryName{
    #use random number as well as localtime to ensure no directory overlap
    my $randomNumber = int( rand(8999) ) + 1000;
    my $directory    = localtime . $randomNumber;
    $directory =~ s/[\W]//g;

    return $directory . '/';
}


sub _createDirectory{
    my $dirName = shift;

    if(defined $dirName){
        #we don't want any permission issues, so don't downgrade the directory permissions
        umask(0);

        #from File::Path
        make_path($dirName) or die "Couldn't create fastaBase $dirName\n";
    }
    else{
        print STDERR "undefined directory name in _createDirectory\n";
        exit(1);
    }

}