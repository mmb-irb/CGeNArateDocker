#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 3){
    print "Usage: perl $0 <curvesAvg> <units> <BasePairParamsAvg.table>\n";
    exit(0);
}

my ($avg_input,$units,$table)=@ARGV;

my @arr = split '/',$avg_input;
my $avg = $arr[$#arr];

my %avgTable = readTable($table);

foreach my $pop (sort keys %avgTable){
	foreach my $bp (sort keys %{$avgTable{$pop}}){
		foreach my $param (sort keys %{$avgTable{$pop}->{$bp}}){
			my $f = $avgTable{$pop}->{$bp}->{$param};
			print "$pop $bp $param: $f\n";
		}
	}
}

my $title = $avg;
$title =~s/.avg.dat//g;
$title = ucfirst($title);
print "Title: -$title-\n";

my $max = -999;
my $min = 999;
my $type = "DNA";
open GNUP, ">$avg_input.gnuplot";
open AVG,"$avg_input";
while(<AVG>){
	chomp;
	# GG   0.01   
	# CC   0.06  

	my ($pair,$value) = split ' ';

	my $bp = substr($pair,0,2);
	$type = "RNA" if ($bp =~ /U/);

	my $dnabsc0 = $avgTable{'Parmbsc0-DNA'}->{$bp}->{$title};
	my $dnacharmm = $avgTable{'Charmm27-DNA'}->{$bp}->{$title};
	my $rnabsc0 = $avgTable{'Parmbsc0-RNA'}->{$bp}->{$title};
	my $rnacharmm = $avgTable{'Charmm27-RNA'}->{$bp}->{$title};

	if($value){
		$max = $value if ($value > $max);
		$min = $value if ($value < $min);
	}
	if($dnabsc0){
		$max = $dnabsc0 if ($dnabsc0 > $max);
		$min = $dnabsc0 if ($dnabsc0 < $min);
	}
	if($rnabsc0){
		$max = $rnabsc0 if ($rnabsc0 > $max);
		$min = $rnabsc0 if ($rnabsc0 < $min);
	}
	if($dnacharmm){
		$max = $dnacharmm if ($dnacharmm > $max);
		$min = $dnacharmm if ($dnacharmm < $min);
	}
	if($rnacharmm){
		$max = $rnacharmm if ($rnacharmm > $max);
		$min = $rnacharmm if ($rnacharmm < $min);
	}

	$dnabsc0 = "?" if (!$dnabsc0);
	$dnacharmm = "?" if (!$dnacharmm);
	$rnabsc0 = "?" if (!$rnabsc0);
	$rnacharmm = "?" if (!$rnacharmm);

	print GNUP "$pair $value $dnabsc0 $dnacharmm $rnabsc0 $rnacharmm\n";

}
close AVG;
close GNUP;

my $scale = ($max - $min) / 5;
print "Max: $max, Min: $min, Scale: $scale\n";
$max += $scale;
$min -= $scale;
print "NewMax: $max, NewMin: $min\n";

my $plot = "gnuplot";
open OUT, ">$plot.tmp";
#print OUT "set xrange[$x_min:$x_max]\n";
print OUT "set yrange[$min:$max]\n";
print OUT "set title 'Base Pair Step Helical Parameters Stiffness Constants: $title'\n";
print OUT "set xlabel 'Sequence Base Pair'\n";
print OUT "set ylabel '$title ($units)'\n";
print OUT "set style line 1 lt 1 lw 2 pt 7 ps 1.5\n";
print OUT "set style line 2 lt 2 lw 1 pt 7 ps 1.5\n";
print OUT "set style line 3 lt 3 lw 1 pt 7 ps 1.5\n";
print OUT "set style line 4 lt 4 lw 1 pt 7 ps 1.5\n";
print OUT "set style line 5 lt 5 lw 1 pt 7 ps 1.5\n";
print OUT "set offsets 0.2, 0.2, $scale, $scale\n";
print OUT "set term png font 'Helvetica,10'\n";
print OUT "set xtics font \"Helvetica,8\"\n";
#print OUT "set term png \n";
print OUT "set output \"$avg_input.png\" \n";
#print OUT "set key out vert center top\n";
#print OUT "plot \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Stiffness Average', \"$avg_input.gnuplot\" using 3:xticlabels(1) with lp ls 2 title '$title Parmbsc0-DNA', \"$avg_input.gnuplot\" using 4:xticlabels(1) with lp ls 3 title '$title Charmm27-DNA', \"$avg_input.gnuplot\" using 5:xticlabels(1) with lp ls 4 title '$title Parmbsc0-RNA', \"$avg_input.gnuplot\" using 6:xticlabels(1) with lp ls 5 title '$title Charmm27-RNA'\n";
if($type eq "DNA"){
print OUT "plot \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Stiffness Average', \"$avg_input.gnuplot\" using 3:xticlabels(1) with lp ls 2 title '$title Parmbsc0-DNA', \"$avg_input.gnuplot\" using 4:xticlabels(1) with lp ls 3 title '$title Charmm27-DNA'\n";
}
else{
print OUT "plot \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Stiffness Average', \"$avg_input.gnuplot\" using 5:xticlabels(1) with lp ls 4 title '$title Parmbsc0-RNA', \"$avg_input.gnuplot\" using 6:xticlabels(1) with lp ls 5 title '$title Charmm27-RNA'\n";
}
print OUT "quit \n";
close OUT;

system ("gnuplot $plot.tmp");
#system ("rm $plot.tmp");


sub readTable{
	my ($file) = @_;

	my $pop;
	my %table;
	open TABLE,"$file";
	while(<TABLE>){
		chomp;
		if(/FF:/){
			# FF: Charmm27-RNA
			my $tag;
			($tag,$pop) = split ':';
			$pop =~s/ //g;
		}
		next if (/^#/);

		#   Twist      Tilt      Roll      Shift     Slide    Rise
		#AA  0.028     0.037     0.020      1.72      2.13     7.64
		#AC  0.036     0.038     0.023      1.28      2.98     8.83

		my ($bp, $twist, $tilt, $roll, $shift, $slide, $rise) = split ' ';

		$table{$pop}->{$bp}->{'Twist'} = $twist;
		$table{$pop}->{$bp}->{'Tilt'} = $tilt;
		$table{$pop}->{$bp}->{'Roll'} = $roll;
		$table{$pop}->{$bp}->{'Shift'} = $shift;
		$table{$pop}->{$bp}->{'Slide'} = $slide;
		$table{$pop}->{$bp}->{'Rise'} = $rise;
	}
	close TABLE;

	return %table;
}
