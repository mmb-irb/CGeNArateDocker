#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 4){
    print "Usage: perl $0 <file.dat> <title> <xlabel> <ylabel>\n";
    exit(0);
}

my ($input,$title,$xlabel,$ylabel)=@ARGV;

my $max = -999;
my $min = 999;
my $len = 0;
open AVG,"$input";
while(<AVG>){
	chomp;
	# K-10 58.6520 0.005
	# S-11 57.6030 0.002

        next if(!$_);

	my ($pair,$value,$stdev) = split ' ';

	if($value ne '?'){
	        my $maxvalue = $value + $stdev;
        	my $minvalue = $value - $stdev;

		$max = $maxvalue if ($maxvalue > $max);
		$min = $minvalue if ($minvalue < $min);
	}
	my $l = length($pair);
	$len = $l if ($l > $len);
}
close AVG;

my $offset = 1;
my $offset_key = 2;
$offset = 2.5 if($len >=5);
my $scale = ($max - $min) / 10;
print "Max: $max, Min: $min, Scale: $scale\n";
$max += $scale;
$min -= $scale;
print "NewMax: $max, NewMin: $min, Key lenght: $len\n";

#$max = $max + 2;
my $plot = "gnuplot";
open OUT, ">$plot.tmp";
#print OUT "set xrange[$x_min:$x_max]\n";
#print OUT "set yrange[$min:$max]\n";
print OUT "set offsets 0.2, 0.2, $scale, $scale\n";
print OUT "set yrange[$min:$max]\n";
print OUT "set title '$title'\n";
print OUT "set xlabel '$xlabel' offset 0, -$offset_key\n";
print OUT "set ylabel '$ylabel'\n";
print OUT "set xtics rotate by -90 offset 0,-$offset \n";
print OUT "set datafile missing \"?\"\n";
#print OUT "set size ratio 0.1\n";
print OUT "set style line 1 lt 1 lw 2 pt 7 ps 1.5\n";
#print OUT "set style line 2 lt 2 lw 2 pt 7 ps 1.5\n";
#print OUT "set term png size 5500, 700\n";
print OUT "set term png font 'Helvetica,10'\n";
print OUT "set xtics font 'Helvetica,8'\n";
print OUT "set output \"$input.png\" \n";
#print OUT "plot \"$input\" using 0:2:3:xticlabels(1) with yerrorlines lp ls 1 title '$title'\n";
print OUT "plot \"$input\" using 0:2:3:xticlabels(1) with yerrorlines notitle, \"$input\" using 2:xticlabels(1) with lp ls 1 title '$title'\n";
print OUT "quit \n";
close OUT;

system ("gnuplot $plot.tmp");
system("cp $plot.tmp $plot.gnuplot.tmp");
#system ("rm $plot.tmp");
