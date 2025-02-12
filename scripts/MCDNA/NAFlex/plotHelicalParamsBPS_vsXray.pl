#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 3){
    print "Usage: perl $0 <curvesAvg> <BasePairStepParamsAvg.table> <curvesXray>\n";
    exit(0);
}

my ($avg_input,$table,$xray)=@ARGV;

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
$title =~s/_avg.dat//g;
$title = ucfirst($title);
print "Title: $title\n";

my $ylabel;
if($avg =~ /shift/ or $avg =~ /slide/ or $avg =~ /rise/){
        $ylabel = "$title (Angstroms)";
}
else{
        $ylabel = "$title (Degrees)";
}

my @xrayValues;
open XRAY,"$xray";
while(<XRAY>){
	chomp;
	# 2     GG   0.01
	my ($r,$pair,$v) = split ' ';
	$xrayValues[$r] = $v;
}
close XRAY;

my $max = -999;
my $min = 999;
my $type = "DNA";
open GNUP, ">$avg_input.gnuplot";
open AVG,"$avg_input";
while(<AVG>){
	chomp;
	# 2	GG   0.01   0.10
	# 3	CC   0.06   0.10

	my ($r,$pair,$value,$stdev) = split ' ';

	$stdev = 0.0 if (!$stdev); # PDB case, just one snapshot.

	my $bp = substr($pair,0,2);
        $type = "RNA" if ($bp =~ /U/);

	my $avgParm;
	my $avgParmMD;
	my $avgParmRNA;

	$avgParm = $avgTable{'All'}->{$bp}->{$title} if(exists $avgTable{'All'}->{$bp}->{$title});
	$avgParmMD = $avgTable{'MD'}->{$bp}->{$title} if(exists $avgTable{'MD'}->{$bp}->{$title});
	$avgParmRNA = $avgTable{'RNA'}->{$bp}->{$title} if(exists $avgTable{'RNA'}->{$bp}->{$title});

	#print "avgTable{'All'}->{$bp}->{$title} = $avgParm\n" ;

        $avgParm = "?" if (! defined $avgParm);
        $avgParmMD = "?" if (! defined $avgParmMD);
        $avgParmRNA = "?" if (! defined $avgParmRNA);

	my $xrayV = $xrayValues[$r];

	print GNUP "$pair $value $stdev $avgParm $avgParmMD $avgParmRNA $xrayV\n";

	my $maxvalue = $value + $stdev;
	my $minvalue = $value - $stdev;

	if($value){
		$max = $maxvalue if ($maxvalue > $max);
		$min = $minvalue if ($minvalue < $min);
	}

	if($xrayV){
		$max = $xrayV if ($xrayV > $max);
		$min = $xrayV if ($xrayV < $min);
	}

	if($type eq 'DNA'){
		if($avgParm && $avgParm ne '?'){
			$max = $avgParm if ($avgParm > $max);
			$min = $avgParm if ($avgParm < $min);
		}
		if($avgParmMD && $avgParmMD ne '?'){
			$max = $avgParmMD if ($avgParmMD > $max);
			$min = $avgParmMD if ($avgParmMD < $min);
		}
	}
	else{
		if($avgParmRNA && $avgParmRNA ne '?'){
			$max = $avgParmRNA if ($avgParmRNA > $max);
			$min = $avgParmRNA if ($avgParmRNA < $min);
		}
	}
}
close AVG;
close GNUP;

my $scale = ($max - $min) / 10;
print "Max: $max, Min: $min, Scale: $scale\n";
$max += $scale;
$min -= $scale;
print "NewMax: $max, NewMin: $min\n";

my $plot = "gnuplot";
open OUT, ">$plot.tmp";
#print OUT "set xrange[$x_min:$x_max]\n";
print OUT "set yrange[$min:$max]\n";
print OUT "set title 'Base Pair Step Helical Parameter: $title'\n";
print OUT "set xlabel 'Sequence Base Pair'\n";
print OUT "set ylabel '$ylabel'\n";
print OUT "set style line 1 lt 1 lw 2 pt 7 ps 1.5\n";
print OUT "set style line 2 lt 2 lw 2 pt 7 ps 1.5\n";
print OUT "set style line 3 lt 3 lw 2 pt 7 ps 1.5\n";
print OUT "set offsets 0.2, 0.2, $scale, $scale\n";
print OUT "set term png font 'Helvetica,10'\n";
print OUT "set xtics font \"Helvetica,8\"\n";
#print OUT "set term png \n";
print OUT "set output \"$avg_input.png\" \n";
#print OUT "set key out vert center top\n";
#print OUT "plot \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Average', \"$avg_input.gnuplot\" using 4:xticlabels(1) with lp ls 2 title '$title ABC Average', \"$avg_input.gnuplot\" using 3:xticlabels(1) with lp ls 3 title '$title X-Ray Average'\n";
if($type eq "DNA"){
print OUT "plot \"$avg_input.gnuplot\" using 0:2:3:xticlabels(1) with yerrorlines ls 1 notitle,  \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Average', \"$avg_input.gnuplot\" using 5:xticlabels(1) with lp ls 2 title '$title ABC Average', \"$avg_input.gnuplot\" using 4:xticlabels(1) with lp ls 3 title '$title X-Ray Average',\"$avg_input.gnuplot\" using 7:xticlabels(1) with lp ls 4 title '$title Experimental X-Ray Structure' \n";
} else{
print OUT "plot \"$avg_input.gnuplot\" using 0:2:3:xticlabels(1) with yerrorlines ls 1 notitle,  \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Average', \"$avg_input.gnuplot\" using 6:xticlabels(1) with lp ls 2 title '$title RNA-MD Average'\n";
}
#print OUT "\"$input\" using 0:$i:$j:xticlabels(1) with yerrorlines ls $cont notitle, \"$input\" using $i:xticlabels(1) with lp ls $cont title '$legends[$cont-1]'";
print OUT "quit \n";
close OUT;

system ("gnuplot $plot.tmp");
#system ("rm $plot.tmp");


sub readTable{
	my ($file) = @_;

	my %table;
	open TABLE,"$file";
	while(<TABLE>){
		chomp;
		next if (/^#/);

		# Population  Step Amount Twist Tilt Roll Shift Slide Rise
		# Naked AA    60  37±8 -1±3  1±4  0.0±0.3 -0.2±0.4 3.3±0.1
		# All   AA   850  34±7 -1±3  3±8  0.1±0.5  0.0±0.6 3.3±0.2
		# MD    AA 25000  35±5 -3±4  0±6 -0.3±0.6 -0.3±0.6 3.3±0.3

		my ($pop, $bp, $amount, $twist, $tilt, $roll, $shift, $slide, $rise) = split ' ';

		my ($twist_v,$twist_stdev) = split '±',$twist;
		my ($tilt_v,$tilt_stdev) = split '±',$tilt;
		my ($roll_v,$roll_stdev) = split '±',$roll;
		my ($shift_v,$shift_stdev) = split '±',$shift;
		my ($slide_v,$slide_stdev) = split '±',$slide;
		my ($rise_v,$rise_stdev) = split '±',$rise;

		$table{$pop}->{$bp}->{'Twist'} = $twist_v;
		$table{$pop}->{$bp}->{'Tilt'} = $tilt_v;
		$table{$pop}->{$bp}->{'Roll'} = $roll_v;
		$table{$pop}->{$bp}->{'Shift'} = $shift_v;
		$table{$pop}->{$bp}->{'Slide'} = $slide_v;
		$table{$pop}->{$bp}->{'Rise'} = $rise_v;

		$table{$pop}->{$bp}->{'Twist_stdev'} = $twist_stdev;
		$table{$pop}->{$bp}->{'Tilt_stdev'} = $tilt_stdev;
		$table{$pop}->{$bp}->{'Roll_stdev'} = $roll_stdev;
		$table{$pop}->{$bp}->{'Shift_stdev'} = $shift_stdev;
		$table{$pop}->{$bp}->{'Slide_stdev'} = $slide_stdev;
		$table{$pop}->{$bp}->{'Rise_stdev'} = $rise_stdev;
	}
	close TABLE;

	return %table;
}
