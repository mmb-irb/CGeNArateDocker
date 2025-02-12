#!/usr/bin/perl
use strict;

# Input Parameters.
my $long=@ARGV;
if ($long != 3){
    print "Usage: perl $0 <NOE.distMat.dat> <NOE.intMat.out.tmp> <NOE.intMat.out>\n";
    exit(0);
}
my ($fnIn, $int, $fnOut) =@ARGV;

open IN,$fnIn;
my $i = 0;
my $line = '';
my @names;
while(<IN>){
	if($i == 0){
		$line = $_;
		$i = 1;
		next;
	}
	my ($v,$rest) = split ' ';
	push(@names,$v);
}
close IN;

open OUT, ">$fnOut";
print OUT $line;
open IN2, "$int";
$i = 0;
while(<IN2>){
	print OUT "$names[$i] $_";
	$i++;
}
close IN2;
close OUT;
