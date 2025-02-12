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
my ($count,$count_hisW,$count_hisC,$bpstep1,$bpstep2);
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
`rm -r ./backbone_torsions/` if (-s "backbone_torsions");
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
}
foreach my $torsion (@backbone_torsions){
	print "Extracting $torsion values ...\n";
	&extract_torsions($torsion);
	# Torsion backbone angle for every base pair (ends removed).
	foreach my $bp (1 .. $seqlen*2){
                open (OUTPUT2, ">backbone_torsions/${torsion}.dat") || die "Cannot open backbone_torsions/${torsion}.dat\n";
                print OUTPUT2 "$bp\t$backbone{$bp}->{$torsion}\n";
                close (OUTPUT2);
        }
}

foreach my $bp (1 .. $seqlen*2){
	print "Calculating epsilon - zeta differences...\n";
        open (OUTPUT2, ">backbone_torsions/diff.dat") || die "Cannot open backbone_torsions/diff.dat\n";
	my $diff = $backbone{$bp}->{epsil}-$backbone{$bp}->{zeta};
	print OUTPUT2 "$bp\t$diff\n";
	close(OUTPUT2);
}

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
	my ($count_canonical_alpha_gamma,$count1,$count_BI,$count2,$count_north,$count_east,$count_south,$count_west,$count3)=(0,0,0,0,0,0,0,0,0);

	my $alpha = $backbone{$bp}->{alpha};
	my $gamma = $backbone{$bp}->{gamma};
	if ($alpha && $gamma){
		$count1++;
		if ((240<$alpha && $alpha<360) && (0<$gamma && $gamma<120)){$count_canonical_alpha_gamma++;}
	}
	my $epsilon = $backbone{$bp}->{epsil};
	my $zeta = $backbone{$bp}->{zeta};
	if ($epsilon && $zeta){
		$count2++;
		my $difference = ($epsilon - $zeta);
		if ($difference < 0){ 
			$count_BI++;
		} 
	}
	my $pucker = $backbone{$bp}->{phase};
	if ($pucker){
                $count3++;
		# Dividing the pseudorotation cycle in four regions.
                if ($pucker < 45 || $pucker > 315){$count_north++;}
                if (45<$pucker && $pucker<135){$count_east++;}
                if (135<$pucker && $pucker<225){$count_south++;}
                if (225<$pucker && $pucker<315){$count_west++;}
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
	open (TORSION,"$file1") || die "Cannot open $file1\n";
        while (<TORSION>){
		my ($bp,$v) = split ' ';
		if ($v < 0){my $a = $v+360;$v = $a;}
		$backbone{$bp}->{$torsion} = $v;
		push (@{$backbone{$torsion}->{$bp}},$v);
	}
	close(TORSION);

	$file2 = `ls *${torsion}C.ser`;
	chomp($file2);
	open (TORSION,"$file2") || die "Cannot open $file2\n";
	while (<TORSION>){
		my ($bp,$v) = split ' ';
		my $bp2=($seqlen*2)-$bp+1;
		if ($v < 0){my $a = $v+360;$v = $a;}
		$backbone{$bp2}->{$torsion} = $v;
		push (@{$backbone{$torsion}->{$bp2}},$v);
	}
	close (TORSION);
}
