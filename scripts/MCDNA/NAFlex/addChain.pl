#!/usr/bin/perl -w
use strict;

my $long=@ARGV;
if ($long != 2){
    print "Usage: perl $0 <inputFile> <outputFile>\n";
    exit(0);
}

my ($file,$out) = @ARGV;

my @ch = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z 0 1 2 3 4 5 6 7 8 9);

open OUT,">$out";

my $i = 0;
my $resn;
my $res;
my $resant= -99999;
my $next = 0;
open FILE,"$file";
while(<FILE>){

	if(/^ATOM/ or /^HETATM/){

		$resn=substr($_,22,5);
		$resn=~s/ //g;
		$res = substr($_,17,4);
		$res =~s/ //g;
		if($resn != $resant and $next) { $i++; $next = 0;}

		my $f = substr($_,0,21);
		my $oldChain = substr($_,21,1);
		my $s = substr($_,22);
		my $chain = $ch[$i];
		if($oldChain eq ' '){
			print OUT "$f$chain$s";
		}
		else{
			print OUT "$_";
		}
	}
	else{
		if (/^TER/){
			my $f = substr($_,0,21);
			print "$f$ch[$i]\n";
		}
		else{
			print OUT "$_";
		}
	}

	if (/^TER/ || $res =~ /3$/) {
		$next = 1;
	}

	$resant = $resn;
}
close FILE;

close OUT;
