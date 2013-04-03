#!/usr/bin/env perl
#for email sending, using the sendEmail (http://caspian.dotconf.net/menu/Software/SendEmail/) program
#Ideally, we would use Email::Sender, but CentOS will not play nice with the requirements for Net::SSLeay

package Interface::Scripts::SendMailWrapper;

use strict;
use warnings;
use IO::File;

#object creation
sub new {
    my ($class) = shift;
    my $self = {};
    bless( $self, $class );
    $self->_initialize(@_);
    return $self;
}


sub _initialize{
    my $self=shift;

    #take in the file for the defaults
    if(@_){
        $self->_setDefaultsFromFile(@_);
    }
    else{
        print STDERR "File of email settings required for Interface::Scripts::SendMailWrapper\n";
        exit(1);
    }
}

#create getters/setters
sub to{
    my $self=shift;
    $self->{'_to'}=shift // return $self->{'_to'};
}

sub from{
    my $self=shift;
    $self->{'_from'}=shift // return $self->{'_from'};
}

sub server{
    my $self=shift;
    $self->{'_server'}=shift // return $self->{'_server'};
}


sub subject{
    my $self=shift;
    $self->{'_subject'}=shift // return $self->{'_subject'};
}

sub message{
    my $self=shift;
    $self->{'_message'}=shift // return $self->{'_message'};
}

sub password{
    my $self=shift;
    $self->{'_password'}=shift // return $self->{'_password'};
}

sub username{
    my $self=shift;
    $self->{'_username'}=shift // return $self->{'_username'};
}

sub footer{
    my $self=shift;
    $self->{'_footer'}=shift // return $self->{'_footer'};
}

sub downloadLink{
    my $self=shift;
    $self->{'_downloadLink'}=shift // return $self->{'_downloadLink'};
}

sub _setDefaultsFromFile{
    my($self)=shift;       
    my $fileName=shift;
    
    my $defaultsFH = IO::File->new('<' . $fileName) or die "cannot open $fileName $!\n";

    while(my $line= $defaultsFH->getline){
            $line =~ s/\R//g;
            my @la = split('\t',$line);
            
            my $parameter = $la[0];
            my $pValue= $la[1];
            
            if($self->can($parameter)){
                $self->$parameter($pValue);
            }                  
    }               
    $defaultsFH->close();   
}


sub sendTheEmail{
    my($self)=shift;
    
    if($self->passesCheck()){
        $self->_createTheMessage();
        my $mailLine = 'sendEmail -t ' . "'" . $self->to() . "'" . ' -f ' . "'" . $self->from() . "'" . ' -s ' . "'" . $self->server() . "'" 
                . ' -xu ' . "'" . $self->username() . "'" . ' -xp ' . "'" . $self->password() . "'" . ' -u ' . "'" . $self->subject() . "'" 
                . ' -m ' . "'" . $self->message() . "'"; 
        
        system("$mailLine");
    }
    else{
            print STDERR "arguments missing for sendEmail!\n";
            exit(1);
    }
}

sub _createTheMessage{
        my($self)=shift;
        
        my $concatString = $self->message() . "\n\n" . $self->downloadLink() . "\n\n" . $self->footer() . "\n";
        $self->message($concatString);
}

sub passesCheck{
        my($self)=shift;
        
        unless(
            $self->to() &&
            $self->from() &&
            $self->server() &&
            $self->username() &&
            $self->password() &&
            $self->subject() &&
            $self->message()
        ){
            print STDERR "Missing required field in Interface::Scripts::SendMailWrapper\n";
            return 0;
        }
        return 1;       
}

1;
