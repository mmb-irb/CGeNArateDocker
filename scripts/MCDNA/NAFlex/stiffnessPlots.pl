#!/usr/bin/perl -w

use strict;

my $outShift = "shift.avg.dat";
my $outSlide = "slide.avg.dat";
my $outRise = "rise.avg.dat";
my $outTilt = "tilt.avg.dat";
my $outRoll = "roll.avg.dat";
my $outTwist = "twist.avg.dat";

my @hps;
$hps[0] = "Shift";
$hps[1] = "Slide";
$hps[2] = "Rise";
$hps[3] = "Tilt";
$hps[4] = "Roll";
$hps[5] = "Twist";

my %bps_avg;
my %bps_cont;
my $cmd = `ls --color=never *.cte`;
foreach my $file (split '\n',$cmd){
	print "$file\n";

	my $bps = uc(substr($file,0,4));

	print "BPS: $bps\n";

	my $cont = 1;
	my @lines;
	open IN,"$file";
	while(<IN>){
		chomp;
		my @array = split '\s+';

		$lines[$cont] = $array[$cont];
print "$array[$cont]\n";
		print "Lines[$cont] = $lines[$cont]\n";
		$cont++;
	}
	close IN;

	for(my $i=0; $i<=5; $i++){
		my $hp = $hps[$i];
		my $j = $i + 1;
		print "Lines[$j]: $lines[$j]\n";
		$bps_avg{$bps}->{$hp} += $lines[$j];
		$bps_cont{$bps}->{$hp}++;
		print "bps_avg {$bps} -> {$hp} += $lines[$j] ($bps_avg{$bps}->{$hp})\n";
	}
}

open SHIFT,">$outShift";
open SLIDE,">$outSlide";
open RISE,">$outRise";
open TILT,">$outTilt";
open ROLL,">$outRoll";
open TWIST,">$outTwist";

foreach my $bps (sort keys %bps_avg){
	foreach my $hp2 (sort keys %{$bps_avg{$bps}}){
		my $val = $bps_avg{$bps}->{$hp2};
		my $cont = $bps_cont{$bps}->{$hp2};
		my $avg = $val;
		$avg = $val / $cont if ($cont);
		print "$val / $cont = $avg\n";
		print SHIFT "$bps $avg\n" if($hp2 eq "Shift");
		print SLIDE "$bps $avg\n" if($hp2 eq "Slide");
		print RISE "$bps $avg\n" if($hp2 eq "Rise");
		print TILT "$bps $avg\n" if($hp2 eq "Tilt");
		print ROLL "$bps $avg\n" if($hp2 eq "Roll");
		print TWIST "$bps $avg\n" if($hp2 eq "Twist");
		print "$bps $hp2 $avg\n";
	}
}

close SHIFT;
close SLIDE;
close RISE;
close TILT;
close ROLL;
close TWIST;
