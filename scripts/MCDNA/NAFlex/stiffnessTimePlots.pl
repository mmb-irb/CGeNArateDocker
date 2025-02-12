#!/usr/bin/perl -w

use strict;

my $cmd = `ls --color=never *.txt`;
foreach my $file (split '\n',$cmd){
	print "$file\n";

	my $bps = uc(substr($file,0,4));

	print "BPS: $bps\n";

	# gatc.4.txt
	my $name = $file;
	$name =~s/\.txt//g;

	my ($bpstep,$num) = split '\.',$name;
	$bpstep = uc($bpstep);
	my $code = "$num-$bpstep";

	mkdir("$code") if(! -e "$code");

	my $outShift = "$code/shift.dat";
	my $outSlide = "$code/slide.dat";
	my $outRise = "$code/rise.dat";
	my $outTilt = "$code/tilt.dat";
	my $outRoll = "$code/roll.dat";
	my $outTwist = "$code/twist.dat";

	open SHIFT,">$outShift";
	open SLIDE,">$outSlide";
	open RISE,">$outRise";
	open TILT,">$outTilt";
	open ROLL,">$outRoll";
	open TWIST,">$outTwist";

	my $cont = 1;
	my @lines;
	open IN,"$file";
	while(<IN>){
		chomp;
		my @array = split '\s+';

		print SHIFT "$cont $array[1]\n";
		print SLIDE "$cont $array[2]\n";
		print RISE "$cont $array[3]\n";
		print TILT "$cont $array[4]\n";
		print ROLL "$cont $array[5]\n";
		print TWIST "$cont $array[6]\n";

		$cont++;
	}
	close IN;

	close SHIFT;
	close SLIDE;
	close RISE;
	close TILT;
	close ROLL;
	close TWIST;

	foreach my $file (`ls --color=never $code/*.dat`){
		chomp($file);
		print "$file\n";

		my $stats = $file;
		$stats =~s/\.dat//g;

		my ($mean,$stdev) = getAvg($file);
		open STDEV,">$stats.stats";
		print STDEV "Mean: $mean, StDev: $stdev\n";
		close STDEV;
	}
}

# getAvg
# Computing Mean and Stdev from a set of values in a 2-rows file.
sub getAvg{

        my ($file) = @_;

        my $mean = 0;
        my $stdev = 0;
        my %info;

        open FILE,"$file";
        while(<FILE>){
                # 1     0.42
                # 2     1.25
                chomp;
                my ($snapshot,$value) = split ' ';
                $info{$snapshot} = $value;
                $mean += $value;
        }
        close FILE;

        my $cont = keys %info;
        $mean /= $cont if($cont);

        foreach my $snp (sort keys %info){
                my $v = $info{$snp};
                my $add = $v - $mean;
                my $add2 = $add * $add;
                $stdev += $add2;
        }
        $stdev /= $cont if($cont);
        $stdev = sqrt($stdev);

        return ($mean,$stdev);
}

