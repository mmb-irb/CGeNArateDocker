#!/usr/bin/perl -w

use strict;

# run_distContactMaps
# Computing Distance Contact Maps from a Nucleic Acids Trajectory.
# Saving output files in $workdir/DistanceContactMaps
sub run_distContactMaps{

	my ($workdir,$prefix,$pdb,$top,$crd) = @_;

	chdir($workdir);
	#mkdir("DistanceContactMaps_$prefix") if(! -e "DistanceContactMaps_$prefix");
	#chdir("DistanceContactMaps_$prefix");
	mkdir("CONTACTS") if(! -e "CONTACTS");
	chdir("CONTACTS");

        ###########  CURVES+  #############

        my ($seq1,$seq2,$firsts,$lasts) = checkPDB($pdb);
        print "SEQ1: $seq1\n";
        print "SEQ2: $seq2\n\n";

	my $revSeq2 = reverse $seq2;
	my $l1 = length($seq1);
	my $l2 = length($seq2);
	my $l = $l1+$l2;

        open SEQR,">seqR.info";
        my $cont = 1;
        foreach my $f (split '',$seq1){ 
                print SEQR "$cont-$f ";
                $cont++;
        }
        foreach my $f (split '',$revSeq2){
                print SEQR "$cont-$f ";
                $cont++;
        }
        print SEQR "\n";
        close SEQR;

	open NUC,">seq.info";
	print NUC "  Strand  1 has  $l1 bases (5'-3'): $seq1\n";
	print NUC "  Strand  2 has  $l2 bases (3'-5'): $seq2\n";
	print NUC "  NucType: DNA\n";
	close NUC;

	my $fullSeq = $seq1.$revSeq2;

        ###########  Computing Distances  #############

	my $output = $prefix."_cgAnalOut";

	print "Getting distances with ptraj program from Ambertools...\n";
	#print "Going to execute: exePtrajCM_Distances($prefix,$output,$pdb,$top,$crd)\n";
	#my $atToTakeIntoAccountInComputingDistances = "C5";
	my $atToTakeIntoAccountInComputingDistances = "C1'";
	exePtrajCM_Distances($prefix,$output,$pdb,$top,$crd,$atToTakeIntoAccountInComputingDistances);

	if(! -s "ptraj.out"){
		print "Cannot find ptraj output file ptraj.out (or is empty). Exiting...\n";
		exit;
	}

	###########  Plotting (R Scripts) - Contact Maps  #############

	print "Plotting Contact Maps...\n\n";
	my $title = "Distances (Angstroms): ";	
	my $out = "distanceMean.contactMap";
	my $ini = 1;
	my $end = $l;
	my $offset = 1;
	`perl $scriptsDir/plotContactMap_R.pl "$title" $ini $end $offset $out `;

	my $outClass = `cat ptraj.out`;

	open OUT, ">nucleicAnalysis.out";
	print OUT "# HBs Output File:\n\n";
	print OUT "$outClass\n";
	close OUT;

        ###########  Building Pdf Files joining png plots together (MuG)  #############

        #`convert *.png ../CONTACTS.pdf`;

	# Cleaning folder...
	`rm *-*.dat`;

	print "DONE!!\n";
}

###########################  SUBROUTINES  ################################

# exePtrajCM_Distances
# Executing ptraj program
sub exePtrajCM_Distances {
	my ($prefix,$out,$pdb,$top,$crd,$at) = @_;

	my %nucleotides;
	my $codeRes;
	my $codeResAnt='';
	my $amber = 1;
	my $newamber = 0;
	my $heteroNucs = 0;
	open PDB,"$pdb";
	while(<PDB>){
		chomp;
		next if !(/^ATOM/);
                my $at = substr($_,12,5);
                $at =~s/ //g;
		if($at=~/^\'/){ 		# New Amber Forcefields (ff10-ff12)
	                $newamber = 1;
			$amber = 0;
		}
		my $res = substr($_,17,3);
		$res =~s/ //g;
		$amber = 0 if($nucleic_codes_charmm{$res});
		my $nres = substr($_,22,5);
		$nres =~s/ //g;
		$codeRes = "$res$nres";
			#print "1Nucleotide: $codeRes - $nres\n";
		if($codeRes ne $codeResAnt){
			$nucleotides{$nres} = $res if ($nucleic_codes{$res});
			$heteroNucs++ if (!$nucleic_codes{$res});
			#print "2Nucleotide: $codeRes - $nres\n";
			$codeResAnt = $codeRes;
		}
	}
	close PDB;

	#ATOM     66  C5' DC  X   3     -13.346   5.182   8.258  0.00  0.00            
	#ATOM     67 1H5' DC  X   3     -13.690   6.195   8.466  0.00  0.00            
	#ATOM     68 2H5' DC  X   3     -14.216   4.542   8.404  0.00  0.00            
	#ATOM     69  C4' DC  X   3     -12.743   5.049   6.792  0.00  0.00            
	#ATOM     70  H4' DC  X   3     -13.592   5.201   6.126  0.00  0.00            
	#ATOM     71  O4' DC  X   3     -11.747   6.110   6.626  0.00  0.00            
	#ATOM     72  C1' DC  X   3     -10.717   5.503   5.933  0.00  0.00            
	#ATOM     73  H1' DC  X   3     -11.100   5.467   4.913  0.00  0.00            
	#ATOM     86  C3' DC  X   3     -12.062   3.724   6.538  0.00  0.00            
	#ATOM     87  H3' DC  X   3     -12.373   3.078   7.360  0.00  0.00            
	#ATOM     88  C2' DC  X   3     -10.590   4.178   6.542  0.00  0.00            
	#ATOM     89 1H2' DC  X   3     -10.197   4.064   7.553  0.00  0.00            
	#ATOM     90 2H2' DC  X   3     -10.069   3.521   5.846  0.00  0.00            
	#ATOM     91  O3' DC  X   3     -12.407   3.125   5.284  0.00  0.00            

	open IN,">ptraj.in";
	print IN "trajin $crd\n";

	my %done;
	foreach my $d (sort {$a <=> $b} keys %nucleotides){
		foreach my $d2 (sort {$a <=> $b} keys %nucleotides){
#	for (my $i=$ini;$i<$end-$offset;$i+=$offset){
#		for (my $j=$i+$offset;$j<=$end;$j+=$offset){
#			next if($done{$d2} or $d == $d2);
			print IN "distance $d-$d2 :$d\@* :$d2\@* out $d-$d2.dat\n"
#			print IN "distance $i-$j :$i\@$at :$j\@$at out $i-$j.dat\n"
		}
#		$done{$d} = 1;
	}
	print IN "analyze statistics all\n";
	close IN;

	print "$amberexe/cpptraj $top < ptraj.in > ptraj.out\n";
	`$amberexe/cpptraj $top < ptraj.in > ptraj.out`;

}

return 1;
