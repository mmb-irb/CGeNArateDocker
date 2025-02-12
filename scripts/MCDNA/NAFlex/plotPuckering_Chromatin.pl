#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 1){
    print "Usage: perl $0 <puckeringAvg>\n";
    exit(0);
}

my ($avg_input)=@ARGV;

my @arr = split '/',$avg_input;
my $avg = $arr[$#arr];

my $title = $avg;
$title =~s/.dat//g;
print "Title: $title\n";

my ($title1,$title2);
if($title =~ /alpha/){
	$title = "Cannonical Alpha-Gamma";
	$title1 = "Gamma";
	$title2 = "Alpha";
}
else{
	$title = "BI/BII Population";
	$title1 = "BII";
	$title2 = "BI";
}

my $totalBasesWC = `wc -l $avg_input`;
my ($totalBases,$f) = split ' ',$totalBasesWC;
$totalBases-- if($avg_input =~ /puckering/);

exit if ($totalBases < 2);

my $numBases = $totalBases / 2;
print "NumBases Strand: $numBases\n";

my $puckering = 0;
my $contBases = 0; 
open GNUP, ">$avg_input.gnuplot";
open AVG,"$avg_input";
while(<AVG>){
	chomp;
	# Puckering:
	#	      north(%) east(%) south(%) west(%)
	# 1	C    0.00   16.83   83.17    0.00

	# Alpha-gamma, BI/BII:
	# 2	G	89.11

	if (/north/){
		$puckering = 1;
		next;
	}

	my @line = split ' ';
	my ($r,$base,$north_or_value,$east,$south,$west) = @line;

        $base = $base."5'" if ($contBases == 0);
        $base = $base."3'" if ($contBases == ($numBases - 1));
        $base = $base."5'" if ($contBases == ($numBases));
        $base = $base."3'" if ($contBases == ($totalBases - 1));
        $base = $base."-".$r;

	if($puckering){
		# Getting cummulative values
		$east += $north_or_value;
		$south += $east;
		$west += $south;
	
		print GNUP "$base $north_or_value $east $south $west\n";
	}
	else{
		print GNUP "$base $north_or_value 100.00\n";
	}

        $contBases++;
        if($contBases == $numBases){
                print GNUP "| 0.00 0.00 0.00 0.00\n";
        }
}
close AVG;
close GNUP;

#my $offset = 3.5;
my $offset = 2.5;
#my $offset_key = 2;
my $offset_key = 2;

#set xlabel 'Nucleotide Sequence' offset 0, -2
#set xtics rotate by 90 offset 0, -2.5 

my $plot = "gnuplot";
open OUT, ">$plot.tmp";
#print OUT "set bmargin 7\n";

print OUT "set terminal unknown\n";

if ($puckering){
	print OUT "set title 'Nucleotide Parameter: Puckering'\n";
	print OUT "set ylabel 'Puckering (%)'\n";
	print OUT "plot \"$avg_input.gnuplot\" using 5:xticlabels(1) with boxes title 'west', \"$avg_input.gnuplot\" using 4:xticlabels(1) with boxes title 'south', \"$avg_input.gnuplot\" using 3:xticlabels(1) with boxes title 'east', \"$avg_input.gnuplot\" using 2:xticlabels(1) with boxes title 'north'\n";
}
else{
	print OUT "set title 'Nucleotide Parameter: $title'\n";
	print OUT "set ylabel '$title (%)'\n";
	print OUT "plot \"$avg_input.gnuplot\" using 3:xticlabels(1) with boxes title '$title1', \"$avg_input.gnuplot\" using 2:xticlabels(1) with boxes title '$title2' \n";
}

print OUT "xspan = GPVAL_DATA_X_MAX - GPVAL_DATA_X_MIN\n";
print OUT "yspan = GPVAL_DATA_Y_MAX - GPVAL_DATA_Y_MIN\n";
print OUT "xequiv = 100\n";
print OUT "yequiv = 400\n";
print OUT "ar = yspan/xspan * xequiv/yequiv\n";
print OUT "ydim = 900\n";
print OUT "xdim = 900/ ar\n";
print OUT "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]\n";
print OUT "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]\n";
#print OUT "set bmargin 5\n";

#print OUT "set xlabel 'Nucleotide Sequence' offset 0, -$offset_key\n";
print OUT "set xlabel 'Nucleotide Sequence' offset 0, -$offset_key\n";
#print OUT "set xtics rotate by -45 offset 0,-$offset \n";
print OUT "set xtics rotate by 90 offset 0, -$offset\n";
print OUT "set key outside\n";

print OUT "set terminal png size xdim,ydim font 'Helvetica,10'\n";
print OUT "set size ratio ar\n";

#print OUT "set term png font 'Helvetica,10'\n";
print OUT "set output \"$avg_input.png\" \n";
print OUT "set style fill solid 0.25 border\n";
print OUT "set boxwidth 0.8\n";
#print OUT "set xtics font \"Helvetica,8\"\n";

print OUT "replot \n";

print OUT "quit \n";
close OUT;

system ("gnuplot $plot.tmp");
#system ("rm $plot.tmp");

