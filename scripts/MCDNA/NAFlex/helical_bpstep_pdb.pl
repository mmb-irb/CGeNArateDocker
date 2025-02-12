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
my ($count,$count_his,$bpstep1,$bpstep2);
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
my $max = $seqlen - 1;
my $max_bps = $max*2;
`rm -r helical_bpstep`;
`mkdir helical_bpstep`;
foreach my $bp (1 .. $max){
#        next if ($bp==$seqlen-1 || $bp==$seqlen);
        $bpstep1 = $fwd[$bp-1].$fwd[$bp];
        $bpstep2 = $rev[$bp].$rev[$bp-1];
        $bpstep{$bp} = $bpstep1.$bpstep2;
}
foreach my $param (keys %bpstephelparam){
	print "Extracting $param values ...\n";
	&extract($param);
	open (AVG,">helical_bpstep/${param}_avg.dat");
	# Base pair step helical parameters for every base pair step (ends removed).
	foreach my $bp (1 .. $max){
		# Base pair step helical parameters.
		print AVG "$bp\t$bpstep{$bp}";
		printf AVG "%7.2f\n",$helical{$param}->{$bp};
	}
	close (AVG);
}

######################## Subroutines ##############################

sub extract {
	my ($param) = @_;
	my $file = `ls *${param}.ser`;
	open PARAM,"$file" || die "Cannot open file\n";
	while (<PARAM>){	
		my ($bp,$v) = split ' ';
		$helical{$param}->{$bp} = $v;
	}
	close (PARAM);
}	
