#!/usr/bin/perl
#Authen-SASL required for gmail
#example taken from http://www.daniweb.com/software-development/perl/threads/230444/error-sending-mail-from-netsmtpssl

use Net::SMTP::SSL;
use strict;

sub send_mail {
my $to = 'sw.laing@gmail.com';
my $subject = 'butterhorn';
my $body = 'hard for the money';
my $from = 'panseq.results@gmail.com';
my $password = '7&phac4$PHAC';
my $smtp;
if (not $smtp = Net::SMTP::SSL->new('smtp.gmail.com',
                            Port => 465,
                            Debug => 1)) {
   die "Could not connect to server\n";
}
$smtp->auth($from, $password)
   || die "Authentication failed!\n";
$smtp->mail($from . "\n");
my @recepients = split(/,/, $to);
foreach my $recp (@recepients) {
    $smtp->to($recp . "\n");
}
$smtp->data();
$smtp->datasend("From: " . $from . "\n");
$smtp->datasend("To: " . $to . "\n");
$smtp->datasend("Subject: " . $subject . "\n");
$smtp->datasend("\n");
$smtp->datasend($body . "\n");
$smtp->dataend();
$smtp->quit;
}
# Send away!
&send_mail();