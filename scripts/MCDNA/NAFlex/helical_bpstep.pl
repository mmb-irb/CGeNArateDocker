#!/usr/bin/perl -w

# Input: Sequence/*ser files from CURVES program.
# Output: base pair step helical parameter for every base pair step and average values for every base pair step helical parameter removing the ends (or a list of bpsteps).
use strict;
my ($fwd,$rev) = @ARGV;
my $long=@ARGV;
if ($long<2){
        print "Please, give me the fwd sequence (5'->3') and the reverse (5'->3')\n";
        exit(0);
}
$fwd = uc($fwd);
$rev = uc($rev);
print "$fwd\t$rev\n";
my $seqlen = length($fwd);
print "Sequence length: $seqlen\n";
my @fwd = split ('', $fwd);
my @rev = split ('', $rev);

my $number = 0;
my $count = 0;
my $count_his = 0;
my ($bpstep1,$bpstep2);
my (%splitfwd,%splitrev,%bpstep,%helical,%helical_his);
foreach my $base (@fwd){
        $number++;
        $splitfwd{"$number"} = $base;
}
foreach my $base (@rev){
        $number++;
        $splitrev{"$number"} = $base;
}
my %bpstephelparam = qw(shift 1 slide 1 rise 1 tilt 1 roll 1 twist 1);
my $max = $seqlen-1;
my $max_bps = $max*2;
`rm -r helical_bpstep`;
`mkdir helical_bpstep`;
foreach my $bp (2 .. $max-1){
        next if ($bp==$seqlen-1 || $bp==$seqlen);
        $bpstep1 = $fwd[$bp-1].$fwd[$bp];
        $bpstep2 = $rev[$bp].$rev[$bp-1];
        $bpstep{$bp} = $bpstep1.$bpstep2;
	`mkdir helical_bpstep/$bp-$bpstep{$bp}`;
        print "$bp\t$bpstep{$bp}\n";
}
foreach my $param (keys %bpstephelparam){
	print "Extracting $param values ...\n";
	&extract($param);
	open (AVG,">helical_bpstep/${param}_avg.dat");
	# Base pair step helical parameters for every base pair step (ends removed).
	foreach my $bp (2 .. $max-1){
		$bpstep1 = $fwd[$bp-1].$fwd[$bp];
		$bpstep2 = $rev[$bp].$rev[$bp-1];
		open (OUTPUT, ">helical_bpstep/$bp-$bpstep{$bp}/${param}.dat");
		foreach my $time (1 .. $count){
			push (@{$helical{$param}->{$bp}},$helical{$bp}->{$param}->{$time});
			# NUCMD Correction, only bps DIFFERENT than shift and tilt with bpstep1 == bpstep2 have params * by -1
			# By now, just commenting it...
			#if ($param eq 'shift' || $param eq 'tilt'){
			#	if ($bpstep1 eq $bpstep2){push (@{$helical{$param}->{$bp}},-$helical{$bp}->{$param}->{$time});}
			#}
			print OUTPUT "$time\t$helical{$bp}->{$param}->{$time}\n";
		}
		close (OUTPUT); 
		open (OUTPUT2, ">helical_bpstep/$bp-$bpstep{$bp}/${param}.his");
		foreach my $time (1 .. $count_his-1){
			#push (@{$helical_his{$param}->{$bp}},$helical_his{$bp}->{$param}->{$time});
			#if ($param eq 'shift' || $param eq 'tilt'){

			#}
			print OUTPUT2 "$helical_his{$bp}->{$param}->{$time}->{'time'}\t$helical_his{$bp}->{$param}->{$time}->{'value'}\n";
		}
		close (OUTPUT2);
		# Average value for the base pair step helical parameters along time.
		my ($avg,$sd)=&stats(@{$helical{$param}->{$bp}});
		print AVG "$bp\t$bpstep{$bp}";
		printf AVG "%7.2f%7.2f\n",$avg,$sd;
	}
	close (AVG);
}
open (SUMM,">helical_bpstep/summary.dat");
foreach my $bp (2 .. $max-1){
	foreach my $param ("shift","slide","rise","tilt","roll","twist"){
		my ($avg,$sd)=&stats(@{$helical{$param}->{$bp}});
		printf SUMM "%7.2f%7.2f",$avg,$sd;
	}
	print SUMM "\n";
}
close (SUMM);
######################## Subroutines ##############################
sub extract {
	my ($param) = @_;
	my $file = `ls *${param}.ser`;
	chomp($file);
	if(-s $file){
		open PARAM,"$file" || die "Cannot open file\n";
		$count=0;
		while (<PARAM>){	
			my @fields = split;
			my $time = $fields[0];
			foreach my $bp (2 .. $max-1){
				next if (!$fields[$bp]);
				$helical{$bp}->{$param}->{$time} = $fields[$bp];
			}
			$count++;
		}
		close (PARAM);
	}

	$file = `ls *${param}.his`;
	chomp($file);
	if(-s $file){
		open PARAM,"$file" || die "Cannot open file\n";
		$count_his=1;
		while (<PARAM>){
			my @fields = split;
			my $time = $fields[0];
			foreach my $bp (2 .. $max){
				next if (!$fields[$bp]);
				$helical_his{$bp}->{$param}->{$count_his}->{'time'} = $fields[0];
				$helical_his{$bp}->{$param}->{$count_his}->{'value'} = $fields[$bp];
			}
			$count_his++;
		}
		close (PARAM);
	}
}	

sub by_number {
    if ($a < $b) {
        -1;
    } elsif ($a == $b) {
        0;
    } elsif ($a > $b) {
        1;
    }
}

sub stats {
	my @in = @_;
        my $number = $#in + 1; #count how many

###total sum calculation###
	
	my $counter = 0;
        my $sd;
        my $total = 0;
        my $mean_rounded;
        while ($counter < $number){
		$total = $total + $in[$counter];
            $counter++;
        }

####mean and standard deviation calc####
	
	my $sum_of_squares = 0;
        if( $number <=1){
                $mean_rounded = $in[0];
                $sd = 0;
        }
        else {
                my $mean = $total / $number;       #calculates average
                my $number2 = $number - 1; #subtract 1 from total for standard deviation calculation
                $counter = 0;
                while ($counter < $number){
                    $sum_of_squares = $sum_of_squares + ($in[$counter] - $mean)**2;
                    $counter++;
                }
                my $standard_deviation = sqrt ($sum_of_squares / $number2);
                $sd = sprintf("%.1f", $standard_deviation);
                $mean_rounded = sprintf("%.2f", $mean); #rounds average to one decimal place with a % on the end
        }
        return ($mean_rounded,$sd);
}

