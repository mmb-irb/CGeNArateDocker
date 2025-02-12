#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 3){
    print "Usage: perl $0 <curvesAvgMD> <curvesAvgExp> <outPlot>\n";
    exit(0);
}

my ($avg_inputMD,$avg_inputExp,$out)=@ARGV;

my @arr = split '/',$avg_inputMD;
my $avg = $arr[$#arr];

my $title = $avg;
$title =~s/_avg.dat//g;
$title = ucfirst($title);
print "Title: $title\n";

my $ylabel;
if($avg =~ /shift/ or $avg =~ /slide/ or $avg =~ /rise/){
        $ylabel = "$title (Angstroms)";
}
elsif($avg =~ /disp/){
	if($avg =~ /xdisp/){
		$title = "X-Displacement";
		$ylabel = "X-Displacement (Angstroms)";
	}
	if($avg =~ /ydisp/){
		$title = "Y-Displacement";
		$ylabel = "Y-Displacement (Angstroms)";
	}
}
elsif($avg =~ /shear/ || $avg =~ /stretch/ || $avg =~ /stagger/){
        $ylabel = "$title (Angstroms)";
}
else{
	$ylabel = "$title (degrees)";
}

my %exp;
open AVG,"$avg_inputExp";
while(<AVG>){
	chomp;
	# 2	GG   0.01   0.10
	my ($r,$pair,$value,$stdev) = split ' ';
	$exp{$r}->{$pair}->{'value'} = $value;	
	$exp{$r}->{$pair}->{'stdev'} = $stdev;	
}
close AVG;

my $max = -999;
my $min = 999;
my $type = "DNA";
open GNUP, ">$avg_inputMD.gnuplot";
open AVG,"$avg_inputMD";
while(<AVG>){
	chomp;
	# 2	GG   0.01   0.10
	# 3	CC   0.06   0.10

	my ($r,$pair,$value,$stdev) = split ' ';

	my $bp = substr($pair,0,2);
        $type = "RNA" if ($bp =~ /U/);

	my $avgExpV = $exp{$r}->{$pair}->{'value'} if(exists $exp{$r}->{$pair}->{'value'});
	my $avgExpStdev = $exp{$r}->{$pair}->{'stdev'} if(exists $exp{$r}->{$pair}->{'stdev'});

        $avgExpV = "?" if (! defined $avgExpV);
        $avgExpStdev = "?" if (! defined $avgExpStdev);

	print GNUP "$pair $value $stdev $avgExpV $avgExpStdev\n";

	my $maxvalue = $value + $stdev;
	my $minvalue = $value - $stdev;

	if($value){
		$max = $maxvalue if ($maxvalue > $max);
		$min = $minvalue if ($minvalue < $min);
	}

	if($avgExpV && $avgExpV ne '?'){
		$max = ($avgExpV + $avgExpStdev) if ( ($avgExpV + $avgExpStdev) > $max);
		$min = ($avgExpV - $avgExpStdev) if ( ($avgExpV - $avgExpStdev) < $min);
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
#print OUT "set style line 3 lt 3 lw 2 pt 7 ps 1.5\n";
print OUT "set offsets 0.2, 0.2, $scale, $scale\n";
print OUT "set term png \n";
print OUT "set output \"$out.dat.png\" \n";

# CGCG 6.08 6.60 7.72
#print OUT "plot \"$avg_inputMD.gnuplot\" using 0:2:3:xticlabels(1) with yerrorlines ls 1 notitle,  \"$avg_inputMD.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Average', \"$avg_inputMD.gnuplot\" using 0:4:5:xticlabels(1) with yerrorlines ls 2 notitle, \"$avg_inputMD.gnuplot\" using 4:xticlabels(1) with lp ls 2 title '$title Experimental'\n";
print OUT "plot \"$avg_inputMD.gnuplot\" using 0:4:5:xticlabels(1) with yerrorlines ls 2 notitle,  \"$avg_inputMD.gnuplot\" using 4:xticlabels(1) with lp ls 2 title 'Experimental $title', \"$avg_inputMD.gnuplot\" using 0:2:3:xticlabels(1) with yerrorlines ls 1 notitle, \"$avg_inputMD.gnuplot\" using 2:xticlabels(1) with lp ls 1 title 'ParmBSC1-MD Average $title'\n";

print OUT "quit \n";
close OUT;

system ("gnuplot $plot.tmp");
