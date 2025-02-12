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
foreach my $bp (1 .. $max-1){
#        next if ($bp==$seqlen-1 || $bp==$seqlen);
        $bpstep1 = $fwd[$bp-1].$fwd[$bp];
        $bpstep2 = $rev[$bp].$rev[$bp-1];
        $bpstep{$bp} = $bpstep1.$bpstep2;
}


foreach my $groo (keys %grooves){
        print "Extracting |$groo| values for ...";
        &extract_grooves($groo);
	open (GROOVES,">grooves/${groo}_avg.dat") || die "Cannot open grooves/${groo}_avg.dat\n";
	# Groove dimensions for every base pair step along time.
	foreach my $bp (1 .. $max-1){
		print ".$bp.";
	        # Values for every groove dimension.
		my $v = $groove{$groo}->{$bp};
                if ($v == 0){
                        print GROOVES "$bp\t$bpstep{$bp}\t?\n";
                }
                else{
                        print GROOVES "$bp\t$bpstep{$bp}\t$v\n";
                }
        }
	print "\n";
	close (GROOVES);
}

######################## Subroutines ##############################
sub extract_grooves {
        my ($param) = @_;
        my $file = `ls *${param}.ser`;
        open PARAM,"$file" || die "Cannot open file\n";
        while (<PARAM>){
		my ($bp,$v) = split ' ';
                $groove{$param}->{$bp} = $v;                
        }
	close (PARAM);
}
