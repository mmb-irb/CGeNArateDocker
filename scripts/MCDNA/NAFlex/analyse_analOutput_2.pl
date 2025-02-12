#! /usr/bin/perl -w

# Input: ficheros generados por el script helical_analysis.pl en /hosts/pluto/ifaustino/al para el calculo de las energias de interaccion entre bases vecinas (energias inter e intracadena).
# Output: stacking energies intra- and inter- for each base of the sequence (central 14-mer) and the last column is the standard deviation. Besides this, computes as well the hbonding energies for every base (hbond.dat)

use strict;

my $working_dir=`pwd`;
chomp($working_dir);
my ($wc1,$wc1_init,$wc1_fin,$wc2,$wc2_init,$wc2_fin) = @ARGV;
my $long=@ARGV;
if ($long<6){
        print "Please, give me the <sequenceWC1> <Initial Res_number> <Final Res_number> <sequenceWC2> <Initial Res_number> <Final Res_number>\n";
        exit(0);
}

$wc1 = uc($wc1);
$wc2 = uc($wc2);
print "$wc1\t$wc2\n";
my $seqlen = length($wc1);
print "Sequence length: $seqlen\n";
my @wc1 = split ('', $wc1);
my @wc2 = split ('', $wc2);
my $number = 0;
my $resnum;
my (%wc1,%wc2);
if ($wc1_init < $wc1_fin) {
	for ($resnum=$wc1_init;$resnum<=$wc1_fin;$resnum++){
        	$wc1{$number} = $resnum;
		$number++;
	}
} else {
	for ($resnum=$wc1_init;$resnum>=$wc1_fin;$resnum--){
        	$wc1{$number} = $resnum;
	        $number++;
      	}
}
$number=0;
if ($wc2_init < $wc2_fin) {
	for ($resnum=$wc2_init;$resnum<=$wc2_fin;$resnum++){
        	$wc2{$number} = $resnum;
		$number++;
	}
} else {
	for ($resnum=$wc2_init;$resnum>=$wc2_fin;$resnum--){
        	$wc2{$number} = $resnum;
	        $number++;
      	}
}

#`rm $working_dir/energies/wc.dat`;
open (TOTAL,">energies/wc2.dat") || die "Cannot open total.dat\n";
print TOTAL "Pair of residues - IntraStrand Stacking - InterStrand Stacking - Total Stacking - H-bonding\n";
# Interaction energies for WC interacting strands. 

foreach my $level (0 .. $seqlen-1){
       my (%stacking,%hbonding);
	next if ($level == 0 || $level == $seqlen-1);
        my $resnum_ourbase = $wc1{$level};
	my $resnum_upperbase = $wc1{$level-1};
#	print "$level\t$resnum_ourbase\t$resnum_upperbase\t$resnum_lowerbase\n";
        my $resnum_comp_base = $wc2{$level};
	my $resnum_upperdiagbase = $wc2{$level-1};
#	print "$level\t$resnum_comp_base\t$resnum_upperdiagbase\t$resnum_lowerdiagbase\n";
	my ($h,$i,$j,$k,$l,$m);

#####   $h   #####   $k   ####
#        |            |
#####   $i   #####   $l   ####

	$i = $resnum_ourbase-1;
	$h = $resnum_upperbase-1;
	$l = $resnum_comp_base-1;
	$k = $resnum_upperdiagbase-1;
        open (OUR, "$working_dir/energies/interact_$i") || die "interact_$i doesn't exist:$!\n";
        while (<OUR>){
                my @ie = split;
                push @{$stacking{$i}->{$h}}, $ie[$h];
                push @{$stacking{$i}->{$k}}, $ie[$k];
		push @{$hbonding{$i}->{$l}}, $ie[$l];
        }
        close (OUR);
        open (COMPLEM,"$working_dir/energies/interact_$l") || die "interact_$l doesn't exist:$!\n";
        while (<COMPLEM>){
                my @ie =split;
                push @{$stacking{$l}->{$k}}, $ie[$k];
                push @{$stacking{$l}->{$h}}, $ie[$h];
        }
        my ($avg_intrastrandA, $sd_intrastrandA) = &stats(@{$stacking{$i}->{$h}});
        my ($avg_intrastrandC, $sd_intrastrandC) = &stats(@{$stacking{$l}->{$k}});
        my ($avg_interstrandE, $sd_interstrandE) = &stats(@{$stacking{$i}->{$k}});
        my ($avg_interstrandG, $sd_interstrandG) = &stats(@{$stacking{$l}->{$h}});
        my ($avg_hbonding, $sd_hbonding) = &stats(@{$hbonding{$i}->{$l}});
	my $total_intrastrand = $avg_intrastrandA + $avg_intrastrandC;
	my $total_interstrand = $avg_interstrandE + $avg_interstrandG;
	my $sd_intrastrand = sqrt(($sd_intrastrandA+$sd_intrastrandC)/4);
	my $sd_interstrand = sqrt(($sd_interstrandE+$sd_interstrandG)/4);
        my $totalstack = $total_intrastrand + $total_interstrand;
        my $sd_totalstack = sqrt(($sd_intrastrand+$sd_interstrand)/2);
	print TOTAL "$wc1{$level} $wc2{$level}\t";
        printf TOTAL "%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",$total_intrastrand,$sd_intrastrand,$total_interstrand,$sd_interstrand,$totalstack,$sd_totalstack,$avg_hbonding,$sd_hbonding;
}

close (TOTAL);

######### Subroutines ########

sub stats {
        my @in = @_;
        my $number = $#in + 1; #count how many
        my $counter = 0;
        my $total = 0;
        while ($counter < $number){
            $total = $total + $in[$counter];
            $counter++;
        }
        my ($mean_rounded, $sd);
        if( $number <=1){
            my $mean_rounded = $in[0], my $sd = 0;
        }
        else {
            my $mean = $total / $number;       #calculates average
            my $number2 = $number - 1; #subtract 1 from total for standard deviation calculation
            $counter = 0;
            my $sum_of_squares = 0;
            while ($counter < $number){
                    $sum_of_squares = $sum_of_squares + ($in[$counter] - $mean)**2;
                    $counter++;
                }
                my $standard_deviation = sqrt ($sum_of_squares / $number2);
                $sd = sprintf("%.1f", $standard_deviation);
                $mean_rounded = sprintf("%.1f", $mean); #rounds average to one decimal place with a % on the end
        }
        return ($mean_rounded,$sd);
}
