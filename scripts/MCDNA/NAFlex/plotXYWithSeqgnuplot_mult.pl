#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long < 6){
    print "Usage: perl $0 <file.dat> <num_lines> <xlabel> <ylabel> <title> <@ legends>\n";
    exit(0);
}

my ($input,$lines,$xlabel,$ylabel,$title,@legends)=@ARGV;

my $max = -999;
my $min = 999;
my $len = 0;
open AVG,"$input";
while(<AVG>){
	chomp;
	# K-10 58.6520
	# S-11 57.6030

	next if (/^#/);
	next if (!$_);

	my ($pair,$value) = split ' ';

	if($value){
		$max = $value if ($value > $max);
		$min = $value if ($value < $min);
	}
	my $l = length($pair);
	$len = $l if ($l > $len);
}
close AVG;

my $offset = 2;
my $offset_key = 2;
$offset = 4.5 if($len >=5);
my $scale = ($max - $min) / 10;
print "Max: $max, Min: $min, Scale: $scale\n";
$max += $scale;
$min -= $scale;
print "NewMax: $max, NewMin: $min, Key lenght: $len\n";

$max = $max + 2;
my $plot = "gnuplot";
open OUT, ">$plot.tmp";
print OUT "set datafile missing \"?\"\n";
#print OUT "set xrange[$x_min:$x_max]\n";
#print OUT "set yrange[$min:$max]\n";
print OUT "set offsets 0.2, 0.2, 2, 0\n";
#print OUT "set yrange[0:$max]\n";
print OUT "set yrange[0:10]\n";
print OUT "set title '$title'\n";
print OUT "set xlabel '$xlabel' offset 0, -$offset_key\n";
print OUT "set ylabel '$ylabel'\n";
print OUT "set xtics rotate by -45 offset 0,-$offset \n";
#print OUT "set size ratio 0.1\n";
foreach (my $i=1; $i<=$lines; $i++){
	print OUT "set style line $i lt $i lw 2 pt 7 ps 1.5\n";
}
#print OUT "set style line 2 lt 2 lw 2 pt 7 ps 1.5\n";
#print OUT "set term png size 5500, 700\n";
print OUT "set term png\n";
print OUT "set output \"$input.png\" \n";
#print OUT "plot \"$input\" using 0:2:3:xticlabels(1) with yerrorlines lp ls 1 title '$title'\n";
#print OUT "plot \"$input\" using 0:2:3:xticlabels(1) with yerrorlines notitle, \"$input\" using 2:xticlabels(1) with lp ls 1 title '$title'";
print OUT "plot ";
my $cont=1;
for (my $i=2;$i<=($lines*2);$i+=2){
my $j = $i + 1;
print OUT "\"$input\" using 0:$i:$j:xticlabels(1) with yerrorlines ls $cont notitle, \"$input\" using $i:xticlabels(1) with lp ls $cont title '$legends[$cont-1]'";
print OUT "," if($j<=($lines*2));
$cont++;
}
print OUT ";\n";
print OUT "quit \n";
close OUT;

system ("gnuplot $plot.tmp");
#system ("rm $plot.tmp");
