#!/usr/bin/perl -w

# Input: Sequence/*ser files from CURVES program.
# Output: base pair helical parameter for every base pair and average values for every base pair helical parameter removing the ends (or a list of bps).
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
my %bpstephelparam = qw(shear 1 stretch 1 stagger 1 buckle 1 propel 1 opening 1);
my $max = $seqlen-1;
`rm -r helical_bp` if (-s "helical_bp");
`mkdir helical_bp` if !( -s "helical_bp");
foreach my $bp (2 .. $max){
        $bpstep{$bp} = $fwd[$bp-1].$rev[$bp-1];
	`mkdir helical_bp/$bp-$bpstep{$bp}`;
        print "$bp\t$bpstep{$bp}\n";
}
foreach my $param (keys %bpstephelparam){
	print "Extracting $param values ...\n";
	&extract($param);
	open (AVG,">helical_bp/${param}_avg.dat");
	# Base pair step helical parameters for every base pair step (ends removed).
	foreach my $bp (2 .. $max){
		open (OUTPUT, ">helical_bp/$bp-$bpstep{$bp}/${param}.dat");
		foreach my $time (1 .. $count){
			push (@{$helical{$param}->{$bp}},$helical{$bp}->{$param}->{$time});
			print OUTPUT "$time\t$helical{$bp}->{$param}->{$time}\n";
		}
		close (OUTPUT); 
		open (OUTPUT2, ">helical_bp/$bp-$bpstep{$bp}/${param}.his");
		foreach my $time (1 .. $count_his-1){
			#push (@{$helical_his{$param}->{$bp}},$helical_his{$bp}->{$param}->{$time});
			print OUTPUT2 "$helical_his{$bp}->{$param}->{$time}->{'time'}\t$helical_his{$bp}->{$param}->{$time}->{'value'}\n";
		}
		close (OUTPUT2); 
		# Average value for the base pair step helical parameters along time.
		my $avg = 0;
		my $sd = 0;
		($avg,$sd)=&stats(@{$helical{$param}->{$bp}});
		print AVG "$bp\t$bpstep{$bp}";
		printf AVG "%7.2f%7.2f\n",$avg,$sd;
	}
	close (AVG);
}

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
			foreach my $bp (2 .. $max){
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

