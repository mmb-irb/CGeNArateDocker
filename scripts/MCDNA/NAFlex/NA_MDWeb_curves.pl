#!/usr/bin/perl -w

use strict;

# run_curves
# Executing Curves+, Canal, and Parsing Scripts
# Saving output files in $workdir/Curves
sub run_curves{

	my ($workdir,$prefix,$pdb,$top,$crd) = @_;

	chdir($workdir);
	mkdir("CURVES")  if(! -e "CURVES");
	chdir("CURVES");
	mkdir("FORCE_CTES")  if(! -e "FORCE_CTES");

	my $naType = nucleicType($pdb);

        print "NA_Type: $naType\n";

	###########  CURVES+  #############

	my $outputCurves = $prefix."_curvesOut";
	my $iniSnap = 0;
	my $endSnap = 0;
	my $jumpSnap = 1;

        my ($seq1,$seq2,$firsts,$lasts) = checkPDB($pdb);
        print "SEQ1: $seq1\n";
        print "SEQ2: $seq2\n\n";

        #my $revSeq2 = reverse $seq2;
        my $l1 = length($seq1);
        my $l2 = length($seq2);

        open NUC,">seq.info";
        print NUC "  Strand  1 has  $l1 bases (5'-3'): $seq1\n";
        print NUC "  Strand  2 has  $l2 bases (3'-5'): $seq2\n";
        print NUC "  NucType: $naType\n";
        close NUC;

        my $fullSeq = $seq1.$seq2;

	foreach my $t (@{$firsts}){
		print "FIRSTS: $t\n";
	}
	foreach my $t (@{$lasts}){
	        #print "LASTS: $t\n";
	}

	my @s1i = split ':',$$firsts[0]; # X-DC5:1
	my @s2i = split ':',$$firsts[1]; 

	my @s1e = split ':',$$lasts[0]; # X-DC3:14
	my @s2e = split ':',$$lasts[1]; 

	my $strand1_init = $s1i[1];
	my $strand1_end = $s1e[1];
	my $strand2_init = $s2i[1];
	my $strand2_end = $s2e[1];

	if(-e "$outputCurves.lis"){
	        print "WARNING: Saving backup of old lis file: $outputCurves.lis to curvesBackup.lis...\n";
		`mv $outputCurves.lis curvesBackup.lis`;
	}
	if(-e "$outputCurves.cda"){
	        print "WARNING: Saving backup of old cda file: $outputCurves.cda to curvesBackup.cda...\n";
	        `mv $outputCurves.cda curvesBackup.cda`;
	}

	print "Executing Curves+ program...\n\n";
	#print "Going to execute: exeCurves($crd,$top,$outputCurves,$iniSnap,$endSnap,$jumpSnap,$strand1_init,$strand1_end,$strand2_init,$strand2_end)\n";
	exeCurves($crd,$top,$outputCurves,$iniSnap,$endSnap,$jumpSnap,$strand1_init,$strand1_end,$strand2_init,$strand2_end);

	if(! -s "$outputCurves.lis"){
		print "Cannot find Curves+ output file $outputCurves.lis (or is empty). Exiting...\n";
		exit;
	}

	###########  CANAL (Histo+Series)  #############

	my $outputCanal = $prefix."_canalOut";

	if(-e "$outputCanal.lis"){
		print "WARNING: Saving backup of old lis file: $outputCanal.lis to canalBackup.lis...\n";
		`mv $outputCanal.lis canalBackup.lis`;
	}

	if(-e "$outputCanal.lis"){
	        print "WARNING: Saving backup of old lis file: $outputCanal.lis to canalBackup.lis...\n";
	        `mv $outputCanal.lis canalBackup.lis`;
	}

	print "Executing Canal program (analysis of Curves+ output)...\n\n";
	exeCanal($outputCanal,$outputCurves,$seq1,1,1,0);

	# Saving canal.in/log
	`mv canal.in canalHistoSeries.in`;
	`mv canal.log canalHistoSeries.log`;

	if(! -s "$outputCanal.lis"){
		print "Cannot find Canal output file $outputCanal.lis (or is empty). Exiting...\n";
		exit;
	}

	###########  CANAL Analysis (Nacho Scripts) - helical_bp  #############

	# Analisis de los intra-base pair parametros.

	print "Analysing Intra-base pair parameters...\n";
	#print "Going to exec helical_bp.pl $seq1 $seq2\n";
	`perl $scriptsDir/helical_bp.pl $seq1 $seq2 &> helical_bp.log`;

	###########  CANAL Analysis (Nacho Scripts) - axis_bp  #############

	# Analisis de los base pair axis parametros.

	print "Analysing Axis base pair parameters...\n";
	#print "Going to exec axis_bp.pl $seq1 $seq2\n";
	`perl $scriptsDir/axis_bp.pl $seq1 $seq2 &> axis_bp.log`;

	###########  CANAL Analysis (Nacho Scripts) - helical_bpstep  #############

	# Analisis de los base pair step parametros.

	print "Analysing Base pair step parameters...\n";
	#print "Going to exec helical_bpstep.pl $seq1 $seq2\n";
	`perl $scriptsDir/helical_bpstep.pl $seq1 $seq2 &> helical_bpstep.log`;

	###########  CANAL Analysis (Nacho Scripts) - backbone torsions  #############

	# Analisis de los dihedros del backbone (alpha, beta, chi, delta, epsilon, gamma, zeta).

	print "Analysing Backbone parameters...\n";
	#print "Going to exec backbone_torsions.pl $seq1 $seq2\n";
	`perl $scriptsDir/backbone_torsions.pl $seq1 $seq2 &> backbone_torsions.log`;

	###########  CANAL Analysis (Nacho Scripts) - grooves  #############

	# Analisis de anchura y profundidad de los surcos.

	print "Analysing Grooves parameters...\n";
	#print "Going to exec grooves.pl $seq1 $seq2\n";
	`perl $scriptsDir/grooves.pl $seq1 $seq2 &> grooves.log`;

	###########  CANAL Analysis (Nacho Scripts) - stiffness  #############

	# Analisis de stiffness.

	print "Analysing Stiffness parameters...\n\n";
	#print "Going to exec force_ctes.pl $seq1 $seq2\n";
	`perl $scriptsDir/force_ctes.pl $seq1 $seq2 &> stiffness.log`;

	###########  Plotting (Gnuplot Scripts) - helical_bp  #############

	print "Building Intra-base pair parameters plots...\n";
	#print "Going to build Gnuplot images\n";
	foreach my $file (`ls --color=never helical_bp/*.dat`){
		chomp($file);
		my $pwd = `pwd`;
		chomp($pwd);
		if (-s $file){
			if ($l1 < 30){
				`perl $scriptsDir/plotHelicalParamsBP.pl $file $HelicalBaseParamsTable &> gnuplot.errors`;
			}
			else{
				`perl $scriptsDir/plotHelicalParamsBP_Chromatin.pl $file $HelicalBaseParamsTable &> gnuplot.errors`;
			}	
		}
	}

	###########  Plotting (Gnuplot Scripts) - axis_bp  #############

	print "Building Axis parameters plots...\n";
	#print "Going to build Gnuplot images\n";
	foreach my $file (`ls --color=never axis_bp/*.dat`){
		chomp($file);
		my $pwd = `pwd`;
		chomp($pwd);
		if (-s $file){
			if ($l1 < 30){
				`perl $scriptsDir/plotHelicalParamsBP.pl $file $HelicalBaseParamsTable &> gnuplot.errors`;
			}
			else{
				`perl $scriptsDir/plotHelicalParamsBP_Chromatin.pl $file $HelicalBaseParamsTable &> gnuplot.errors`;
			}	
		}
	}

	###########  Plotting (Gnuplot Scripts) - helical_bpstep  #############

	print "Building Base-pair step parameters plots...\n";
	foreach my $file (`ls --color=never helical_bpstep/*.dat`){
		chomp($file);
		next if ($file =~ /summary/);

		my $filePDB = "../PDB/CURVES/$file";
		
		if (-s $file && -s $filePDB){
			`perl $scriptsDir/plotHelicalParamsBPS_vsXray.pl $file $HelicalParamsTable $filePDB&> gnuplot.errors`;
			#print "Computing: perl $scriptsDir/plotHelicalParamsBPS.pl $file $HelicalParamsTable\n";
		}
		elsif (-s $file){
			if ($l1 < 30){
				`perl $scriptsDir/plotHelicalParamsBPS.pl $file $HelicalParamsTable &> gnuplot.errors`;
			#print "Computing: perl $scriptsDir/plotHelicalParamsBPS.pl $file $HelicalParamsTable\n";
			}
			else{
				`perl $scriptsDir/plotHelicalParamsBPS_Chromatin.pl $file $HelicalParamsTable &> gnuplot.errors`;
			}	
		}
	}

	###########  Plotting (Gnuplot Scripts) - backbone_torsions  #############

	print "Building Backbone parameters plots...\n";
	foreach my $file (`ls --color=never backbone_torsions/*.dat`){
		chomp($file);
		next if ($file =~ /summary/);
		if (-s $file){
			if ($l1 < 30){
				`perl $scriptsDir/plotPuckering.pl $file &> gnuplot.errors`;
			}
			else{
				`perl $scriptsDir/plotPuckering_Chromatin.pl $file &> gnuplot.errors`;
			}
		}
	}
	if ( (-s "${outputCanal}_alphaW.ser") && (-s "${outputCanal}_gammaW.ser") ){
		if ($l1 < 30){
			`perl $scriptsDir/plotAlphaGammaGGplot.pl ${outputCanal}_alphaW.ser ${outputCanal}_gammaW.ser AlphaGamma.dat 'Alpha-Gamma Plot' 'Alpha' 'Gamma' &> plotR_AlphaGamma.log`;
		}
		else{
			`perl $scriptsDir/plotAlphaGammaGGplot_Chromatin.pl ${outputCanal}_alphaW.ser ${outputCanal}_gammaW.ser AlphaGamma.dat 'Alpha-Gamma Plot' 'Alpha' 'Gamma' &> plotR_AlphaGamma.log`;
		}
		`mv AlphaGamma.dat* backbone_torsions/`;
	}
	if ( (-s "${outputCanal}_epsilW.ser") && (-s "${outputCanal}_zetaW.ser") ){
		if ($l1 < 30){
			`perl $scriptsDir/plotAlphaGammaGGplot.pl ${outputCanal}_epsilW.ser ${outputCanal}_zetaW.ser EpsilonZeta.dat 'Epsilon-Zeta Plot' 'Epsilon' 'Zeta' &> plotR_EpsilonZeta.log`;
		}
		else{
			`perl $scriptsDir/plotAlphaGammaGGplot_Chromatin.pl ${outputCanal}_epsilW.ser ${outputCanal}_zetaW.ser EpsilonZeta.dat 'Epsilon-Zeta Plot' 'Epsilon' 'Zeta' &> plotR_EpsilonZeta.log`;
		}
		`mv EpsilonZeta.dat* backbone_torsions/`;
	}


	###########  Plotting (Gnuplot Scripts) - grooves  #############

	print "Building Grooves parameters plots...\n";
	foreach my $file (`ls --color=never grooves/*.dat`){
		chomp($file);
		if (-s $file){
			#`perl $scriptsDir/plotHelicalParamsBP.pl $file $HelicalBaseParamsTable &> gnuplot.errors`;
			if ($l1 < 30){
				`perl $scriptsDir/plotHelicalParamsBP.pl $file $HelicalBaseParamsTable &> gnuplot.errors`;
			#print "Computing: perl $scriptsDir/plotHelicalParamsBPS.pl $file $HelicalParamsTable\n";
			}
			else{
				`perl $scriptsDir/plotHelicalParamsBP_Chromatin.pl $file $HelicalBaseParamsTable &> gnuplot.errors`;
			}	
		}
	}

	###########  Building Packed Curves/Canal Serial & Histogram Output files  #############

	`gzip -f *.cda`;
	my $logTarSer = `tar -zcvf $outputCanal.ser.tgz *.ser`;
	my $logTarHis = `tar -zcvf $outputCanal.his.tgz *.his`;

	###########  Building Output File  #############

	my $curvesOut = `cat $outputCurves.lis`;
	my $canalOut = `cat $outputCanal.lis`;
	open OUT, ">nucleicAnalysis.out";
	print OUT "# Curves+ Output File:\n\n";
	print OUT "$curvesOut\n";
	print OUT "# Canal Output File:\n\n";
	print OUT "$canalOut\n";
	close OUT;

	print "DONE!!\n";
}

return 1;
