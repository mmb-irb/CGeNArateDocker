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
	my %residues = exePtrajCM_Distances($prefix,$output,$pdb,$top,$crd,$atToTakeIntoAccountInComputingDistances);

	if(! -s "ptraj.out"){
		print "Cannot find ptraj output file ptraj.out (or is empty). Exiting...\n";
		exit;
	}


	foreach my $chain (sort {$a <=> $b} keys %residues){

		open SEQR,">seqR.$chain.info";

		my $cont = 1;
		foreach my $residue (sort {$a <=> $b} keys %{$residues{$chain}}){
       
			my $f = $residues{$chain}->{$residue};
        	        #print SEQR "$residue-$f ";	# If we want to see the real residue number from the PDB file
        	        print SEQR "$cont-$f ";		# If we want to see the residue number from 1 to n
                	$cont++;
		}
        	print SEQR "\n";
        	close SEQR;
        }

	###########  Plotting (R Scripts) - Contact Maps  #############

	print "Plotting Contact Maps...\n\n";

	my $title = "Distances (Angstroms): ";	
	my $out = "distanceMean.contactMap";

	# NUC-NUC
	
	my $ini = 1;
	my $end = $l;
	my $offset = 1;
	chdir("NUC-NUC");
	#`cp ../seqR.info .`;
	#`perl $scriptsDir/plotContactMap_R.pl "$title" $ini $end $offset $out `;
	`python $scriptsDir/plotContactMap_Python.py $end`;

	# END TO END distances: distances from nucleotide number 1 to the last one.
	#my $end2end = `ls --color=never -rlht 1-$l1.dat | tail -1 | awk '{print \$NF}'`;
	#chomp($end2end);
	#`cp $end2end end2end.dat`;
	`cp 1-$l1.dat end2end.dat`;

	# LOG file
	my $outClass = `cat ptraj.out`;

        # Cleaning folder...
        #`rm *-*.dat`;
	`find . -type f -name '*-*.dat' -exec rm -f {} +`;

	# PROT-PROT

	chdir("../PROT-PROT");

#my $aa = $pp{'A'};
#my $bb = $pp{'B'};
#
#my ($p1,$p2) = max_min(%{$aa});
#my ($p3,$p4) = max_min(%{$bb});

	my %done;
	foreach my $chain (sort {$a <=> $b} keys %residues){

		my $resCh = $residues{$chain};
		my ($iniProt,$endProt) = max_min(%{$resCh});

		foreach my $chain2 (sort {$a <=> $b} keys %residues){

			next if ($chain eq $chain2);
			next if ($done{"$chain2-$chain"});
			
			my $resCh2 = $residues{$chain2};
			my ($iniProt2,$endProt2) = max_min(%{$resCh2});

			#print "perl $scriptsDir/plotContactMap_prots_R.pl \"$title\" $iniProt $endProt $iniProt2 $endProt2 $offset prot$chain.$chain2.$out \n";
			#`perl $scriptsDir/plotContactMap_prots_R.pl "$title" $iniProt $endProt $iniProt2 $endProt2 $offset prot$chain.$chain2.$out `;
			print "python $scriptsDir/plotContactMap_prots_Python.py \"$title\" $iniProt $endProt $iniProt2 $endProt2 $offset prot$chain-$chain2.$out \n";
			`python $scriptsDir/plotContactMap_prots_Python.py "$title" $iniProt $endProt $iniProt2 $endProt2 $offset prot$chain-$chain2.$out`;

			$done{"$chain-$chain2"} = 1;
		}
	}
	# LOG file
	$outClass .= `cat ptraj.out`;

        # Cleaning folder...
        #`rm *-*.dat`;
	unlink for grep { /^\d+-\d+\.dat$/ } glob("*-*.dat");

	# PROT-NUC

	chdir("../PROT-NUC");

	foreach my $chain2 (sort {$a <=> $b} keys %residues){

		my $resCh = $residues{$chain2};
		my ($ini2,$end2) = max_min(%{$resCh});

		#print "perl $scriptsDir/plotContactMap_prots_R.pl \"$title\" $ini $end $ini2 $end2 $offset nuc.$chain2.$out \n";
		#`perl $scriptsDir/plotContactMap_prots_R.pl "$title" $ini $end $ini2 $end2 $offset nuc.$chain2.$out `;
		#print "python $scriptsDir/plotContactMap_Python.py $end\n";
		#`python $scriptsDir/plotContactMap_Python.py $end`;
		print "python $scriptsDir/plotContactMap_prots_Python.py \"$title\" $ini $end $ini2 $end2 $offset nuc-$chain2.$out \n";
		`python $scriptsDir/plotContactMap_prots_Python.py "$title" $ini $end $ini2 $end2 $offset nuc-$chain2.$out`;
	}

	# LOG file
	$outClass .= `cat ptraj.out`;

        # Cleaning folder...
        #`rm *-*.dat`;
	unlink for grep { /^\d+-\d+\.dat$/ } glob("*-*.dat");

	chdir("..");

	# LOG file
	open OUT, ">nucleicAnalysis.out";
	print OUT "# Contact Distances Output File:\n\n";
	print OUT "$outClass\n";
	close OUT;

        ###########  Building Pdf Files joining png plots together (MuG)  #############

        #`convert *.png ../CONTACTS.pdf`;

	print "DONE!!\n";
}

###########################  SUBROUTINES  ################################

# exePtrajCM_Distances
# Executing ptraj program
sub exePtrajCM_Distances {
	my ($prefix,$out,$pdb,$top,$crd,$at) = @_;

	my %nucleotides;
	my %residues;
	my $codeRes;
	my $codeResAnt='';
	my $amber = 1;
	my $newamber = 0;
	#my $heteroNucs = 0;
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
		my $ch = substr($_,21,1);
		$ch =~s/ //g;
		$codeRes = "$res$nres";
		#print "$_\n";
		#print "1Nucleotide: $codeRes - $nres - $ch\n";
		if($codeRes ne $codeResAnt){
			$nucleotides{$nres} = $res if ($nucleic_codes{$res});
			$residues{$ch}->{$nres} = $res if (!$nucleic_codes{$res});
			#$heteroNucs++ if (!$nucleic_codes{$res});
			#print "2Nucleotide: $ch - $codeRes - $nres\n";
		
			#print "Nucleotide: nucleotides{$nres} = $nucleotides{$nres}\n" if ($nucleic_codes{$res});
			#print "Residue: residues{$ch}->{$nres} =$residues{$ch}->{$nres}\n" if (!$nucleic_codes{$res});

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

	# NUC - NUC distances

	mkdir("NUC-NUC") if (! -s "NUC-NUC");

	open IN,">ptraj.in";
	print IN "trajin $crd\n";
	
	my %done;
	foreach my $d (sort {$a <=> $b} keys %nucleotides){
		foreach my $d2 (sort {$a <=> $b} keys %nucleotides){
			#print IN "distance $d-$d2 :$d\@* :$d2\@* out NUC-NUC/$d-$d2.dat\n"
			#print IN "distance $d-$d2 :$d\@C5 :$d2\@C5 out NUC-NUC/$d-$d2.dat\n"
			print IN "distance $d-$d2 :$d\@C1' :$d2\@C1' out NUC-NUC/$d-$d2.dat\n" if (!$done{"$d-$d2"} and !$done{"$d2-$d"});
			#print IN "distance $d-$d2 ::A:$d\@C1' ::B:$d2\@C1' out NUC-NUC/$d-$d2.dat\n"
			$done{"$d-$d2"} = 1;
			$done{"$d2-$d"} = 1;
		}
	}
	print IN "analyze statistics all\n";
	close IN;

	#`$amberexe/cpptraj $top < ptraj.in > ptraj.out`;
	#`$amberexe/cpptraj14 $pdb < ptraj.in > ptraj.out`;
	`cpptraj $pdb < ptraj.in > ptraj.out`;

	# PROT - PROT distances

	mkdir("PROT-PROT") if (! -s "PROT-PROT");

	open IN,">ptraj.prots.in";
	print IN "trajin $crd\n";

	my %done;
	foreach my $chain (sort {$a <=> $b} keys %residues){
		foreach my $chain2 (sort {$a <=> $b} keys %residues){
			next if ($chain eq $chain2);
			foreach my $d (sort {$a <=> $b} keys %{$residues{$chain}}){
				foreach my $d2 (sort {$a <=> $b} keys %{$residues{$chain2}}){
					print IN "distance $d-$d2 :$d\@CA :$d2\@CA out PROT-PROT/$d-$d2.dat\n" if (!$done{"$d-$d2"});
					#print IN "distance $d-$d2 ::$chain:$d\@CA ::$chain2:$d2\@CA out PROT-PROT/$d-$d2.dat\n" if (!$done{"$d-$d2"});
					#print IN "distance $d-$d2 :$d\@NC :$d2\@NC out PROT-PROT/$d-$d2.dat\n" if (!$done{"$d-$d2"});
					#print IN "distance $d-$d2 :$d\@* :$d2\@* out PROT-PROT/$d-$d2.dat\n" if (!$done{"$d-$d2"});
					$done{"$d-$d2"} = 1;
				}
			}
		}
	}
	print IN "analyze statistics all\n";
	close IN;

	#`$amberexe/cpptraj $top < ptraj.prots.in > ptraj.prots.out`;
	#`$amberexe/cpptraj14 $pdb < ptraj.prots.in > ptraj.prots.out`;
	`cpptraj $pdb < ptraj.prots.in > ptraj.prots.out`;

	# PROT (pilota) - PROT (pilota) distances

	mkdir("PROTpi-PROTpi") if (! -s "PROTpi-PROTpi");

	open IN,">ptraj.prots_pi.in";
	print IN "trajin $crd\n";

	my %done;
	foreach my $chain (sort {$a <=> $b} keys %residues){
		foreach my $chain2 (sort {$a <=> $b} keys %residues){
			next if ($chain eq $chain2);
					print IN "distance $chain-$chain2 ::$chain  ::$chain2 out PROTpi-PROTpi/$chain-$chain2.dat\n" if (!$done{"$chain-$chain2"});
					$done{"$chain-$chain2"} = 1;
		}
	}
	print IN "analyze statistics all\n";
	close IN;

	#`$amberexe/cpptraj $top < ptraj.prots.in > ptraj.prots.out`;
	#`$amberexe/cpptraj14 $pdb < ptraj.prots.in > ptraj.prots.out`;
	`cpptraj $pdb < ptraj.prots_pi.in > ptraj.prots_pi.out`;

	# PROT - NUC distances

	mkdir("PROT-NUC") if (! -s "PROT-NUC");

	open IN,">ptraj.prot-nuc.in";
	print IN "trajin $crd\n";

	undef %done;
	foreach my $d (sort {$a <=> $b} keys %nucleotides){
		foreach my $chain (sort {$a <=> $b} keys %residues){
			foreach my $d2 (sort {$a <=> $b} keys %{$residues{$chain}}){
				print IN "distance $d-$d2 :$d\@* :$d2\@CA out PROT-NUC/$d-$d2.dat\n" if (!$done{"$d-$d2"});
				#print IN "distance $d-$d2 $d\@* ::$chain:$d2\@CA out PROT-NUC/$d-$d2.dat\n" if (!$done{"$d-$d2"});
				#print IN "distance $d-$d2 :$d\@* :$d2\@NC out PROT-NUC/$d-$d2.dat\n" if (!$done{"$d-$d2"});
				#print IN "distance $d-$d2 :$d\@* :$d2\@* out PROT-NUC/$d-$d2.dat\n" if (!$done{"$d-$d2"});
				$done{"$d-$d2"} = 1;
			}
		}
	}
	print IN "analyze statistics all\n";
	close IN;

	#`$amberexe/cpptraj $top < ptraj.prot-nuc.in > ptraj.prot-nuc.out`;
	#`$amberexe/cpptraj14 $pdb < ptraj.prot-nuc.in > ptraj.prot-nuc.out`;
	`cpptraj $pdb < ptraj.prot-nuc.in > ptraj.prot-nuc.out`;

	return %residues;
}

sub max_min {

	my %res = @_;

	my $min = 999999;
	my $max = -999999;
	
	foreach my $v (sort {$a <=> $b} keys %res){
		$min = $v if ($v < $min);
		$max = $v if ($v > $max);
	}

	return ($min,$max);
}

return 1;
