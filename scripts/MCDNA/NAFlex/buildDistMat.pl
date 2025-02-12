#!/usr/bin/perl
use strict;

# Input Parameters.
my $long=@ARGV;
if ($long != 2){
    print "Usage: perl $0 <NOE.stats> <outputDistanceMatrix>\n";
    exit(0);
}
my ($fnIn, $fnOut) =@ARGV;

my %hash;
my %allKeys;
open IN,$fnIn;
while(<IN>){
	chomp;
	# 1 H2pp-H8_1-2 4.5057 0.9130
	# 1 H3p-H8_1-2 4.8438 0.3973
	next if (/^#/);
	my ($snp,$code,$dist,$stdev) = split ' ';
	my ($p1,$p2,$p3) = split '-',$code;
	$p2 =~s/_.*//g;
	my $proton1 = "$snp-$p1";
	my $proton2 = "$p3-$p2";

	$hash{$proton1}->{$proton2} = $dist;
	$hash{$proton2}->{$proton1} = $dist;
	$allKeys{$proton1} = 1;
	$allKeys{$proton2} = 1;
#	print "$proton1 $proton2 $dist\n";
}
close IN;

#foreach my $i (sort {$a <=> $b} keys %hash){
#	foreach my $j (sort {$a <=> $b} keys %{$hash{$i}}){
#		my $v = $hash{$i}->{$j};
#		print "$i $j $v\n";
#	}
#}

open OUT, ">$fnOut";
open OUT2, ">$fnOut.forInt";
foreach my $i (sort {$a <=> $b} keys %allKeys){
	print OUT "$i ";
}
print OUT "\n";

foreach my $i (sort {$a <=> $b} keys %allKeys){
	print OUT "$i ";
	foreach my $j (sort {$a <=> $b} keys %allKeys){

#my ($num1,$name1) = split '-',$i;
#my ($num2,$name2) = split '-',$j;

#		if($num1 eq $num2){
#			print OUT "0.0 ";
#			next;
#		}

		if(exists($hash{$i}->{$j})){
			my $v = $hash{$i}->{$j};
			print OUT "$v ";
			print OUT2 "$v ";
		}
		else{
#			if($i eq $j){
				print OUT "0.0 ";
				print OUT2 "0.0 ";
#			}
#			else{
#				print OUT "100.0 ";
#			}
		}
	}
	print OUT "\n";
	print OUT2 "\n";
}
close OUT;
close OUT2;
