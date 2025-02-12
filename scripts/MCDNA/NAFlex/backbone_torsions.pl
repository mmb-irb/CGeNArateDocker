#!/usr/bin/perl -w

# Input: Sequence/*ser files from CURVES program.
# Output: base pair step dihedral angles evolution along time for every base pair step and percentages for canonical alpha/gamma and BI conformation for every base pair step.

# V1.1: dihedral angles converted to 0º-360º scale (Curves+ does it but from -180ª to +180º scale).
# V2.2: dihedral angles rearranged like in Curves+, that is, by residue and not by step.

use strict;
my ($fwd,$rev) = @ARGV;
my $long=@ARGV;
if ($long<2){
        print "Please, give me the fwd sequence (5'->3') and the reverse (3'->5')\n";
        exit(0);
}
$fwd = uc($fwd);
$rev = uc($rev);
#$rev = reverse($rev);
print "$fwd\t$rev\n";
my $seqlen = length($fwd);
print "Sequence length: $seqlen\n";
my @fwd = split ('', $fwd);
my @rev = split ('', $rev);
my $sequence = $fwd.$rev;
my @sequence = split ('',$sequence);

my $number = 0;
my $count = 0;
my $count_hisW = 0;
my $count_hisC = 0;
my ($bpstep1,$bpstep2);
my (%splitfwd,%splitrev,%bpstep,%backbone,%backbone_his,%dihedral);
foreach my $base (@fwd){
        $number++;
        $splitfwd{"$number"} = $base;
}
foreach my $base (@rev){
        $number++;
        $splitrev{"$number"} = $base;
}
my @backbone_torsions = qw(alpha beta chi gamma epsil phase zeta); 
`rm -r ./backbone_torsions/`;
`mkdir backbone_torsions`;
my $max = $seqlen-1;
my $max_bps = $max*2;
foreach my $bp (1 .. $seqlen*2){
#	next if ($bp==$seqlen-1 || $bp==$seqlen);
	if ($bp < $seqlen+1) {
#		$bpstep1 = $fwd[$bp-1].$fwd[$bp];
#		$bpstep2 = $rev[$bp].$rev[$bp-1];
#		$bpstep{$bp} = $bpstep1.$bpstep2;
		$bpstep{$bp} = $fwd[$bp-1];
	}
        else {
                my $bp2 = $bp-$seqlen; 
#		$bpstep1 = $rev[$max-$bp2+1].$rev[$max-$bp2];
#               $bpstep2 = $fwd[$max-$bp2].$fwd[$max-$bp2+1];
#		$bpstep{$bp} = $bpstep1.$bpstep2;
		$bpstep{$bp} = $rev[$max-$bp2+1];
        }
	print "$bp\t$bpstep{$bp}\n";
	`mkdir backbone_torsions/$bp-$bpstep{$bp}`;
}
foreach my $torsion (@backbone_torsions){
	print "Extracting $torsion values ...\n";
	&extract_torsions($torsion);
	# Torsion backbone angle for every base pair (ends removed).
	foreach my $bp (1 .. $seqlen*2){
                open (OUTPUT2, ">backbone_torsions/$bp-$bpstep{$bp}/${torsion}.dat") || die "Cannot open backbone_torsions/$bp-$bpstep{$bp}/${torsion}.dat\n";
                foreach my $time (1 .. $count){
                        print OUTPUT2 "$time\t$backbone{$bp}->{$torsion}->{$time}\n";
                }
                close (OUTPUT2);
        }

	foreach my $bp (1 .. $seqlen*2){
                open (OUTPUT2, ">backbone_torsions/$bp-$bpstep{$bp}/${torsion}.his") || die "Cannot open backbone_torsions/$bp-$bpstep{$bp}/${torsion}.his\n";
		my $count_his = $count_hisW;
		$count_his = $count_hisC if ($bp > $seqlen);
                foreach my $time (1 .. $count_his-1){
                        print OUTPUT2 "$backbone_his{$bp}->{$torsion}->{$time}->{'time'}\t$backbone_his{$bp}->{$torsion}->{$time}->{'value'}\n";
                }
                close (OUTPUT2);
        }
}

foreach my $bp (1 .. $seqlen*2){
	print "Calculating epsilon - zeta differences...\n";
        open (OUTPUT2, ">backbone_torsions/$bp-$bpstep{$bp}/diff.dat") || die "Cannot open backbone_torsions/$bp-$bpstep{$bp}/diff.dat\n";
        foreach my $time (1 .. $count){
		my $diff = $backbone{$bp}->{epsil}->{$time}-$backbone{$bp}->{zeta}->{$time};
		print OUTPUT2 "$time\t$diff\n";
	}
	close(OUTPUT2);
}

open SUMM,">backbone_torsions/summary.dat";
foreach my $bp (1..$seqlen*2){
	foreach my $param ("alpha","beta","chi","epsil","gamma","phase","zeta"){
		my ($avg,$sd) = &statsAngles(@{$backbone{$param}->{$bp}});
		printf SUMM "%7.2f%7.2f",$avg,$sd;
	}	
	print SUMM "\n";
}
close SUMM;

#############
# Calculating population for canonical and alternative alpha_gamma distributions and BI/BII distributions. 
#############

print "Calculating canonical alpha/gamma populations ...\n";
print "... and BI -> BII transitions ...\n";
print "... and puckering populations for every base pair step.\n";
open (ALPHA,">backbone_torsions/canonical_alpha_gamma.dat") || die "Cannot open canonical file\n";
open (TRANS,">backbone_torsions/BI_population.dat") || die "Cannot open BI population file\n";
open (PUCKER,">backbone_torsions/puckering.dat") || die "Cannot open puckering file\n";
print PUCKER "\t      north(%) east(%) south(%) west(%)\n";
foreach my $bp (1 .. $seqlen*2){
#	next if ($bp == $max || $bp == $max+1);
	#my ($count_canonical_alpha_gamma,$count1,$count_BI,$count2,$count_north,$count_east,$count_south,$count_west,$count3);
	my ($count_canonical_alpha_gamma,$count1,$count_BI,$count2,$count_north,$count_east,$count_south,$count_west,$count3) = (0,0,0,0,0,0,0,0,0);
	foreach my $time (1 .. $count){
		my $alpha = $backbone{$bp}->{alpha}->{$time};
		my $gamma = $backbone{$bp}->{gamma}->{$time};
		if ($alpha && $gamma){
			$count1++;
			if ((240<$alpha && $alpha<360) && (0<$gamma && $gamma<120)){$count_canonical_alpha_gamma++;}
#			if ((120<$alpha && $alpha<240) && (120<$gamma && $gamma<240)){$count_alternative_alpha_gamma++;}
		}
		my $epsilon = $backbone{$bp}->{epsil}->{$time};
		my $zeta = $backbone{$bp}->{zeta}->{$time};
		if ($epsilon && $zeta){
			$count2++;
			my $difference = ($epsilon - $zeta);
			if ($difference < 0){ 
				$count_BI++;
			} 
		}
		my $pucker = $backbone{$bp}->{phase}->{$time};
		if ($pucker){
                       $count3++;
			# Dividing the pseudorotation cycle in four regions.
                       if ($pucker < 45 || $pucker > 315){$count_north++;}
                       if (45<$pucker && $pucker<135){$count_east++;}
                       if (135<$pucker && $pucker<225){$count_south++;}
                       if (225<$pucker && $pucker<315){$count_west++;}
                }
	}
	my $canonical_ag_population = $count_canonical_alpha_gamma*100/$count1;
        my $canonical_ag_population_rounded = sprintf ("%.2f",$canonical_ag_population);
        my $BI_population = $count_BI*100/$count2;
        my $BI_population_rounded = sprintf ("%.2f",$BI_population);
        my $north_population = $count_north*100/$count3;
        my $north_population_rounded = sprintf ("%.2f",$north_population);
        my $east_population = $count_east*100/$count3;
        my $east_population_rounded = sprintf ("%.2f",$east_population);
        my $south_population = $count_south*100/$count3;
        my $south_population_rounded = sprintf ("%.2f",$south_population);
        my $west_population = $count_west*100/$count3;
        my $west_population_rounded = sprintf ("%.2f",$west_population);
	print ALPHA "$bp\t$bpstep{$bp}\t$canonical_ag_population_rounded\n";
	print TRANS "$bp\t$bpstep{$bp}\t$BI_population_rounded\n";
	print PUCKER "$bp\t$bpstep{$bp}";
	printf PUCKER "%8.2f%8.2f%8.2f%8.2f\n",$north_population_rounded,$east_population_rounded,$south_population_rounded,$west_population_rounded;
}
close (ALPHA);
close (TRANS);
close (PUCKER);

#######################
#     Subroutines     #
#######################

sub extract_torsions {
	my $torsion=shift;
	my ($file1,$file2);
	$file1=`ls *${torsion}W.ser`;
        chomp($file1);
	if (-s $file1){
		open (TORSION,"$file1") || die "Cannot open $file1\n";
		$count=0;
	        while (<TORSION>){
        	        my @fields = split;
                	my $time = $fields[0];
	                foreach my $bp (1 .. $seqlen){
				next if (!$fields[$bp]);
				if ($fields[$bp] < 0){my $a = $fields[$bp]+360;$fields[$bp] = $a;}
				$backbone{$bp}->{$torsion}->{$time} = $fields[$bp];
				push (@{$backbone{$torsion}->{$bp}},$fields[$bp]);
	                 }
			$count++;
		}
		close(TORSION);
	}

	$file2 = `ls *${torsion}C.ser`;
	chomp($file2);
	if (-s $file2){
		open (TORSION,"$file2") || die "Cannot open $file2\n";
		while (<TORSION>){
	                my @fields = split;
        	        my $time = $fields[0];
                	foreach my $bp (1 .. $seqlen){
				next if (!$fields[$bp]);
				my $bp2=($seqlen*2)-$bp+1;
				if ($fields[$bp] < 0){my $a = $fields[$bp]+360;$fields[$bp] = $a;}
				$backbone{$bp2}->{$torsion}->{$time} = $fields[$bp];
				push (@{$backbone{$torsion}->{$bp2}},$fields[$bp]);
	        	}
		}
		close (TORSION);
	}

        $file1 = `ls *${torsion}W.his`;
	chomp($file1);
	if (-s $file1){
	        open TORSION,"$file1" || die "Cannot open $file1\n";
        	$count_hisW=1;
	        while (<TORSION>){
        	        my @fields = split;
                	my $time = $fields[0];
	                foreach my $bp (1 .. $seqlen){
        	                next if (!$fields[$bp]);
                	        $backbone_his{$bp}->{$torsion}->{$count_hisW}->{'time'} = $fields[0];
	                        $backbone_his{$bp}->{$torsion}->{$count_hisW}->{'value'} = $fields[$bp];
        	        }
                	$count_hisW++;
	        }
        	close (TORSION);
	}

        $file2 = `ls *${torsion}C.his`;
	chomp($file2);
	if (-s $file2){
	        open TORSION,"$file2" || die "Cannot open $file2\n";
	        $count_hisC=1;
        	while (<TORSION>){
	                my @fields = split;
        	        my $time = $fields[0];
                	foreach my $bp (1 .. $seqlen){
	                        next if (!$fields[$bp]);
				my $bp2=($seqlen*2)-$bp+1;
                	        $backbone_his{$bp2}->{$torsion}->{$count_hisC}->{'time'} = $fields[0];
	                        $backbone_his{$bp2}->{$torsion}->{$count_hisC}->{'value'} = $fields[$bp];
        	        }
                	$count_hisC++;
		}
		close (TORSION);
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

sub statsAngles {
	my @in = @_;
        my $number = $#in + 1; #count how many

	my $pi = 3.14159265;

	my @sinAngle;
	my @cosAngle;
	for (my $i=0; $i <= $#in; $i++) {

		# Paso de polar/grados a cartesiano/radianes
		#awk '{print sin($1*3.14159265/180), cos($1*3.14159265/180)}' $name.dat > $name.temp.cart.rad.dat
		my $angle = $in[$i];
		my $s = sin($angle * $pi / 180);
		my $c = cos($angle * $pi / 180);
		$sinAngle[$i] = $s;
		$cosAngle[$i] = $c;
	}

	# Calculo el promedio
	# awk -f averageMOD.awk $name.temp.cart.rad.dat > $name.temp.cart.rad.avg.dat
	my ($meanSin,$sdSin) = &stats(@sinAngle);
	my ($meanCos,$sdCos) = &stats(@cosAngle);

	# Paso de cartesiano/radianes a polar/grados
	# awk '{print (atan2($1,$2)*180/3.14159265)}' $name.temp.cart.rad.avg.dat > $name.avg.dat
	my $mean = atan2($meanSin,$meanCos) * 180 / $pi;
	my $sd = atan2($sdSin,$sdCos) * 180 / $pi;

	$mean = sprintf("%.2f", $mean); 
	$sd = sprintf("%.1f", $sd);

	return ($mean,$sd);
}
