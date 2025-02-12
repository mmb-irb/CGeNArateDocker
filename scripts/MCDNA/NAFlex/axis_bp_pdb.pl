#!/usr/bin/perl -w

# Input: Sequence/*ser files from CURVES program.
# Output: base pair axis parameter for every base pair and average values for every base pair axis parameter removing the ends (or a list of bps).
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
my %bpstephelparam = qw(xdisp 1 ydisp 1 tip 1 inclin 1);
my $max = $seqlen;
`rm -r axis_bp` if(-s "axis_bp");
`mkdir axis_bp`;
foreach my $bp (1 .. $max){
        $bpstep{$bp} = $fwd[$bp-1].$rev[$bp-1];
}
foreach my $param (keys %bpstephelparam){
	print "Extracting $param values ...\n";
	&extract($param);
	open (AVG,">axis_bp/${param}_avg.dat");
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
