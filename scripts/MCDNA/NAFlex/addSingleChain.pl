#!/usr/bin/perl -w
use strict;

my $long=@ARGV;
if ($long != 2){
    print "Usage: perl $0 <inputFile> <outputFile>\n";
    exit(0);
}

my ($file,$out) = @ARGV;

my $single_ch = 'A';

open OUT,">$out";

my $i = 0;
my $cont = 0;
my $resnAnt = 0;
my %done;
open FILE,"$file";
while(<FILE>){

	# Skipping Ions
	next if (/Cl-/ or /Na+/ or /K+/);

	if(/^ATOM/ or /^HETATM/){

		my $resn=substr($_,22,5);
		$resn=~s/ //g;
		$cont++ if ($resn != $resnAnt);
		my $mon=substr($_,17,3);
		$mon=~s/ //g;
		
		if ($mon =~ /5$/ and !$done{$cont} and $cont > 1){
			$i++;
			$done{$cont} = 1;
		}

		my $f = substr($_,0,21);
		my $oldChain = substr($_,21,1);
		my $s = substr($_,22);
		#my $chain = $ch[$i];
		print OUT "$f$single_ch$s";
		#if($oldChain eq ' '){
		#	print OUT "$f$chain$s";
		#}
		#else{
		#	print OUT "$_";
		#}
		$resnAnt = $resn;
	}
	else{
		print OUT "$_" if (!/^TER/);
	}
	
	if(/^TER/ and !$done{$cont}){
		$cont=0;
		$i++;
	}
}
print OUT "TER\n";
close FILE;

close OUT;
