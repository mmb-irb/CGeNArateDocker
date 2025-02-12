#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 2){
    print "Usage: perl $0 <curvesAvg> <BasePairParamsAvg.table>\n";
    exit(0);
}

my ($avg_input,$table)=@ARGV;

my @arr = split '/',$avg_input;
my $avg = $arr[$#arr];

my %avgTable = readTable($table);

foreach my $type (sort keys %avgTable){
        foreach my $bp (sort keys %{$avgTable{$type}}){
	        foreach my $param (sort keys %{$avgTable{$type}->{$bp}}){
			my $f = $avgTable{$type}->{$bp}->{$param};
			print "$type $bp $param: $f\n";
	        }
	}
}

my $title = $avg;
$title =~s/_avg.dat//g;
$title = ucfirst($title);
print "Title: $title\n";

my $numBases;
my $ylabel;
if($avg =~ /min|maj/){ 
	$numBases = "Tetramer";
	if($avg =~ /majw/){
		$title = "Major Groove";
		$ylabel = "Major Groove Width (Angstroms)";
	}
	if($avg =~ /majd/){
		$title = "Major Groove";
		$ylabel = "Major Groove Depth (Angstroms)";
	}
	if($avg =~ /minw/){
		$title = "Minor Groove";
		$ylabel = "Minor Groove Width (Angstroms)";
	}
	if($avg =~ /mind/){
		$title = "Minor Groove";
		$ylabel = "Minor Groove Depth (Angstroms)";
	}
}
elsif($avg =~ /disp/){
	$numBases = "Base";
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
        $numBases = "Pair";
        $ylabel = "$title (Angstroms)";
}
else{
	$numBases = "Pair";
	$ylabel = "$title (degrees)";
}

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

	$type = "RNA" if ($pair =~ /U/);

        my $avgParm;
        my $avgParmRNA;

        $avgParm = $avgTable{'DNA'}->{$pair}->{$title} if(exists $avgTable{'DNA'}->{$pair}->{$title});
        $avgParmRNA = $avgTable{'RNA'}->{$pair}->{$title} if(exists $avgTable{'RNA'}->{$pair}->{$title});

        #print "avgTable{$pair}->{$title} = $avgParm\n" ;

        $avgParm = "?" if (! defined $avgParm);
        $avgParmRNA = "?" if (! defined $avgParmRNA);

	print GNUP "$pair $value $stdev $avgParm $avgParmRNA\n";

	my $maxvalue = -99;
	my $minvalue = 99;
	$maxvalue = $value + $stdev if ($stdev);
	$minvalue = $value - $stdev if ($stdev);

	if($value){
		$max = $maxvalue if ($maxvalue > $max);
		$min = $minvalue if ($minvalue < $min);
	}

	if($type eq 'DNA'){
        	if($avgParm && $avgParm ne '?'){
                	$max = $avgParm if ($avgParm > $max);
	                $min = $avgParm if ($avgParm < $min);
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
print OUT "set title 'Base $numBases Helical Parameter: $title'\n";
print OUT "set xlabel 'Sequence Base $numBases'\n";
print OUT "set ylabel '$ylabel'\n";
print OUT "set style line 1 lt 1 lw 2 pt 7 ps 1.5\n";
print OUT "set offsets 0.2, 0.2, $scale, $scale\n";
print OUT "set style line 2 lt 2 lw 2 pt 7 ps 1.5\n";
print OUT "set term png font 'Helvetica,10'\n";
print OUT "set xtics font \"Helvetica,8\"\n";
#print OUT "set term png \n";
print OUT "set output \"$avg_input.png\" \n";
#print OUT "plot \"$avg.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title MD Average', \"$avg.gnuplot\" using 3:xticlabels(1) with lp ls 2 title '$title ABC Average'\n";
#print OUT "plot \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title MD Average'\n";

if($type eq "DNA"){
	print OUT "plot \"$avg_input.gnuplot\" using 0:2:3:xticlabels(1) with yerrorlines ls 1 notitle, \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Average', \"$avg_input.gnuplot\" using 4:xticlabels(1) with lp ls 2 title '$title ABC Average'\n";
}
else{
	print OUT "plot \"$avg_input.gnuplot\" using 0:2:3:xticlabels(1) with yerrorlines ls 1 notitle, \"$avg_input.gnuplot\" using 2:xticlabels(1) with lp ls 1 title '$title User-MD Average', \"$avg_input.gnuplot\" using 5:xticlabels(1) with lp ls 2 title '$title RNA-MD Average'\n";
}
print OUT "quit \n";
close OUT;

system ("gnuplot $plot.tmp");
#system ("rm $plot.tmp");

sub readTable{
        my ($file) = @_;

        my %table;
	my $type = 'DNA';
        open TABLE,"$file";
        while(<TABLE>){
                chomp;

		$type = 'RNA' if (/RNA/);

                next if (/^#/);

		# BasePair Shear Stretch Stagger Buckle Propeller Opening Xdisp Ydisp Inclination Tip 
		# AT  0.10±0.29  0.03±0.12  0.07±0.42  3.7±11.7 -12.9±9.0  3.1±5.5  -1.45±0.82  0.10±0.53  6.7±5.5  0.2±5.0

                my ($bp, $shear, $stretch, $stagger, $buckle, $propeller, $opening, $xdisp, $ydisp, $inclination, $tip) = split ' ';

                my ($shear_v,$shear_stdev) = split '±',$shear;
                my ($stretch_v,$stretch_stdev) = split '±',$stretch;
                my ($stagger_v,$stagger_stdev) = split '±',$stagger;
                my ($buckle_v,$buckle_stdev) = split '±',$buckle;
                my ($propeller_v,$propeller_stdev) = split '±',$propeller;
                my ($opening_v,$opening_stdev) = split '±',$opening;
                my ($xdisp_v,$xdisp_stdev) = split '±',$xdisp;
                my ($ydisp_v,$ydisp_stdev) = split '±',$ydisp;
                my ($inclination_v,$inclination_stdev) = split '±',$inclination;
                my ($tip_v,$tip_stdev) = split '±',$tip;

                $table{$type}->{$bp}->{'Shear'} = $shear_v;
                $table{$type}->{$bp}->{'Stretch'} = $stretch_v;
                $table{$type}->{$bp}->{'Stagger'} = $stagger_v;
                $table{$type}->{$bp}->{'Buckle'} = $buckle_v;
                $table{$type}->{$bp}->{'Propel'} = $propeller_v;
                $table{$type}->{$bp}->{'Opening'} = $opening_v;
                $table{$type}->{$bp}->{'Xdisp'} = $xdisp_v;
                $table{$type}->{$bp}->{'Ydisp'} = $ydisp_v;
                $table{$type}->{$bp}->{'Inclin'} = $inclination_v;
                $table{$type}->{$bp}->{'Tip'} = $tip_v;

                $table{$type}->{$bp}->{'Shear_stdev'} = $shear_stdev;
                $table{$type}->{$bp}->{'Stretch_stdev'} = $stretch_stdev;
                $table{$type}->{$bp}->{'Stagger_stdev'} = $stagger_stdev;
                $table{$type}->{$bp}->{'Buckle_stdev'} = $buckle_stdev;
                $table{$type}->{$bp}->{'Propel_stdev'} = $propeller_stdev;
                $table{$type}->{$bp}->{'Opening_stdev'} = $opening_stdev;
                $table{$type}->{$bp}->{'Xdisp_stdev'} = $xdisp_stdev;
                $table{$type}->{$bp}->{'Ydisp_stdev'} = $ydisp_stdev;
                $table{$type}->{$bp}->{'Inclin_stdev'} = $inclination_stdev;
                $table{$type}->{$bp}->{'Tip_stdev'} = $tip_stdev;

        }
        close TABLE;

        return %table;
}
		
