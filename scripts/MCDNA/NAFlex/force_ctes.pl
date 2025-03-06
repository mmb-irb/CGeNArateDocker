#!/usr/bin/perl
use strict;
my ($fwd,$rev) = @ARGV;
my $long=@ARGV;
if ($long<2){
        print "Please, give me the fwd sequence (5'->3') and the reverse (5'->3')\n";
        exit(0);
}
    
#my $softForces = "/home/NAFlex/soft/force_ctes";
my $softForces = "/app/Scripts/MCDNA/NAFlex/soft/force_ctes";

&generate_files_forces($fwd,$rev);

my @files = qw(shift slide rise tilt roll twist);


sub generate_files_forces {
my ($fwd,$rev) = @_;
$fwd = lc($fwd);
$rev = lc($rev);
print "$fwd\t$rev\n";
my $seqlen = length($fwd);
print "Sequence length: $seqlen\n";
my @fwd = split ('', $fwd);
my @rev = split ('', $rev);

my $number = 0;
my ($count,$bpstep1,$bpstep2);
my (%splitfwd,%splitrev,%bpstep,%backbone,%dihedral);
foreach my $base (@fwd){
        $number++;
        $splitfwd{"$number"} = $base;
}
foreach my $base (@rev){
        $number++;
        $splitrev{"$number"} = $base;
}

my $max = $seqlen-1;
my $max_bps = $max*2;
foreach my $bp (2 .. $max-1){
        next if ($bp==$seqlen-1 || $bp==$seqlen);
        $bpstep1 = $fwd[$bp-1].$fwd[$bp];
        $bpstep2 = $rev[$bp].$rev[$bp-1];
        $bpstep{$bp} = $bpstep1.$bpstep2;
        print "$bp\t$bpstep{$bp}\n";
}
   
   my (@temp,$ok,$line);
   for my $i (2 .. $max-1){
	next if ($i==$max || $i==$max+1);	# Remove from analysis ending base pairs.
	my $shift=`ls *shift.ser`;
	my $slide=`ls *slide.ser`;
	my $rise=`ls *rise.ser`;
	my $tilt=`ls *tilt.ser`;
	my $roll=`ls *roll.ser`;
	my $twist=`ls *twist.ser`;
     open(SH,"<$shift") || die "Cannot open shift file\n";
     open(SL,"<$slide");
     open(RI,"<$rise");
     open(TI,"<$tilt");
     open(RO,"<$roll");
     open(TW,"<$twist");
     $bpstep1 = $fwd[$i-1].$fwd[$i];
     $bpstep2 = $rev[$i].$rev[$i-1];
     open(OUT,">FORCE_CTES/$bpstep{$i}.$i.txt");
#     open(OUT,">tetranucleotides/FORCE_CTES/$bp$complem{$bp}.$i.txt");

	$ok = 1; # MDWeb: We are only interested in the Sequence Base Pairs, not in the reverse ones.
#     $ok = 0;
#     $ok = 1 if ($bpstep1 eq $bpstep2);

     open(OUT2,">FORCE_CTES/$bpstep2$bpstep1.$i.txt") if (!$ok);
#     open(OUT2,">tetranucleotides/FORCE_CTES/$complem{$bp}$bp.$i.txt") if (!$ok);

	# Helical Parameters Order for the next loop: Shift,Slide,Rise,Tilt,Roll,Twist
     while (<SH>) {
       $line = "";
       @temp = split;
       printf OUT "%7.2f",$temp[$i];

       #$line .= sprintf("%7.2f",$temp[$i]*(-1));
       # Adam/Fede (21-12-2016)
       # Removed change of sign (*-1) because it was already considered in computing raw values.
       # When Nacho wrote the script, he needed to fix this here, now is fixed before getting here.
       # SHIFT
       $line .= sprintf("%7.2f",$temp[$i]);
       $_ = <SL>;
       @temp = split;
       printf OUT "%7.2f",$temp[$i];
       $line .= sprintf("%7.2f",$temp[$i]);
       $_ = <RI>;
       @temp = split;
       printf OUT "%7.2f",$temp[$i];
       $line .= sprintf("%7.2f",$temp[$i]);
       $_ = <TI>;
       @temp = split;
       printf OUT "%7.2f",$temp[$i];

       #$line .= sprintf("%7.2f",$temp[$i]*(-1));
       # Adam/Fede (21-12-2016)
       # Remo.ved change of sign (*-1) because it was already considered in computing raw values.
       # When Nacho wrote the script, he needed to fix this here, now is fixed before getting here.
       # TILT
       $line .= sprintf("%7.2f",$temp[$i]);
       $_ = <RO>;
       @temp = split;
       printf OUT "%7.2f",$temp[$i];
       $line .= sprintf("%7.2f",$temp[$i]);
       $_ = <TW>;
       @temp = split;
       printf OUT "%7.2f\n",$temp[$i];
       $line .= sprintf("%7.2f\n",$temp[$i]);
       print OUT2 $line if (!$ok);
       print OUT $line if ($ok);
     }
     close OUT;
     close (SH);
     close (SL);
     close (RI);
     close (TI);
     close (RO);
     close (TW);
     close (OUT2) if (!$ok);
     &calculate_forces("$bpstep{$i}.$i");
     &calculate_forces("$bpstep2$bpstep1.$i") if (!$ok);
   }
       
}

sub calculate_forces {
    my $file = shift;

    print "$file\n";
    open (TMP,">temp.csh");
    print TMP<<EOF;
$softForces/forces.adam -x FORCE_CTES/$file.txt -n 6 -helical > FORCE_CTES/$file.cte
mv average.txt FORCE_CTES/$file.av
EOF
    close(TMP);
    #system("csh temp.csh;rm -f temp.csh");
    system("sh temp.csh;rm -f temp.csh");
}
