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
my (%splitfwd,%splitrev,%bpstep,%groove,%groove_his);
foreach my $base (@fwd){
        $number++;
        $splitfwd{"$number"} = $base;
}
foreach my $base (@rev){
        $number++;
        $splitrev{"$number"} = $base;
}
my %grooves = qw(minw 1 majw 1 mind 1 majd 1);
my $max = $seqlen-1;
my $max_bps = $max*2;
`rm -r grooves`;
`mkdir grooves`;
foreach my $bp (2 .. $max-1){
        next if ($bp==$seqlen-1 || $bp==$seqlen);
        $bpstep1 = $fwd[$bp-1].$fwd[$bp];
        $bpstep2 = $rev[$bp].$rev[$bp-1];
        $bpstep{$bp} = $bpstep1.$bpstep2;
	`mkdir grooves/$bp-$bpstep{$bp}`;
        print "$bp\t$bpstep{$bp}\n";
}


foreach my $groo (keys %grooves){
        print "Extracting |$groo| values for ...";
        &extract_grooves($groo);
	open (GROOVES,">grooves/${groo}_avg.dat") || die "Cannot open grooves/${groo}_avg.dat\n";
	# Groove dimensions for every base pair step along time.
	foreach my $bp (2 .. $max-1){
		print ".$bp.";
		open (OUTPUT, ">grooves/$bp-$bpstep{$bp}/${groo}.dat") || die "Cannot open grooves/$bp-$bpstep{$bp}/${groo}.dat\n";
                foreach my $time (1 .. $count){
                        push (@{$groove{$groo}->{$bp}},$groove{$bp}->{$groo}->{$time}) if ($groove{$bp}->{$groo}->{$time} ne "NaN");
        		print OUTPUT "$time\t$groove{$bp}->{$groo}->{$time}\n";
                }
                close (OUTPUT);
		open (OUTPUT2, ">grooves/$bp-$bpstep{$bp}/${groo}.his");
		foreach my $time (1 .. $count_his-1){
			#push (@{$helical_his{$param}->{$bp}},$helical_his{$bp}->{$param}->{$time});
			print OUTPUT2 "$groove_his{$bp}->{$groo}->{$time}->{'time'}\t$groove_his{$bp}->{$groo}->{$time}->{'value'}\n";
		}
		close (OUTPUT2);
	        # Average value for every groove dimension along time.
        	my ($avg,$sd)=&stats(@{$groove{$groo}->{$bp}});
                if ($avg == 0){
                        print GROOVES "$bp\t$bpstep{$bp}\t?\t?\n";
                }
                else{
                        print GROOVES "$bp\t$bpstep{$bp}\t$avg\t$sd\n";
                }
        }
	print "\n";
	close (GROOVES);
}

######################## Subroutines ##############################
sub extract_grooves {
        my ($param) = @_;
        my $file = `ls *${param}.ser`;
	chomp($file);
	if (-s $file){
	        open PARAM,"$file" || die "Cannot open file\n";
		$count=0;
	        while (<PARAM>){
        	        my @fields = split;
                	my $time = $fields[0];
	                foreach my $bp (2 .. $max-1){
        	        	$groove{$bp}->{$param}->{$time} = $fields[$bp];
                	}
			$count++;
        	}
		close (PARAM);
	}

	$file = `ls *${param}.his`;
	chomp($file);
	if (-s $file){
		open PARAM,"$file" || die "Cannot open file\n";
		$count_his=1;
		while (<PARAM>){
			my @fields = split;
			my $time = $fields[0];
			foreach my $bp (2 .. $max){
				next if (!$fields[$bp]);
				$groove_his{$bp}->{$param}->{$count_his}->{'time'} = $fields[0];
				$groove_his{$bp}->{$param}->{$count_his}->{'value'} = $fields[$bp];
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

