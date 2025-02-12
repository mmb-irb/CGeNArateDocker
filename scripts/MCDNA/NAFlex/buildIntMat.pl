#!/usr/bin/perl
use strict;

# Input Parameters.
my $long=@ARGV;
if ($long != 2){
    print "Usage: perl $0 <NOE.int> <outputIntMatrix>\n";
    exit(0);
}
my ($fnIn, $fnOut) =@ARGV;

my $firstLine = 1;
my @legend;
open OUT, ">$fnOut";
open IN,$fnIn;
while(<IN>){
	chomp;
	# H2pp-H8_1-2 4.5057 0.9130
	# H3p-H8_1-2 4.8438 0.3973
	next if (/^#/);
	if($firstLine){
		@legend = split ' ';
		$firstLine = 0;
		next;
	}
	my ($code,@values) = split ' ';
	my ($num,$prot) = split '-',$code;

	for (my $i=0; $i<=$#values; $i++){
		my $v = $values[$i];
		my $k = $legend[$i];
		my ($numK,$protK) = split '-',$k;
		if($numK == $num){
			print OUT "0.0 ";
		}
		else{
			print OUT "$v ";
		}
	}
	print OUT "\n";
}
close IN;
close OUT;
