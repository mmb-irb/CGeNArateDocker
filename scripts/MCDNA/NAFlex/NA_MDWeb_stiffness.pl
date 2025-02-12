#!/usr/bin/perl -w

use strict;

# run_stiffness
# Executing Curves+, Canal, and Parsing Scripts
# Saving output files in $workdir/Curves
sub run_stiffness{

	my ($workdir,$prefix,$pdb,$top,$crd) = @_;

	my ($seq1,$seq2,$firsts,$lasts);

	my $outputCurves = $prefix."_curvesOut";
	my $outputCanal = $prefix."_canalOut";

	chdir($workdir);
	#mkdir("Stiffness_$prefix")  if(! -e "Stiffness_$prefix");
	#chdir("Stiffness_$prefix");
	mkdir("STIFFNESS")  if(! -e "STIFFNESS");
	chdir("STIFFNESS");

	# Checking if Curves Analysis has already been computed.
	print "Checking if Curves Analysis has already been computed...\n";
	my $analysis = `ls --color=never $workdir`;
	#if ($analysis =~ /(Curves_$prefix)/){
	if ($analysis =~ /(CURVES)/){
		my $dirCurves = $1;

		print "Curves Analysis has already been computed, using Curves files to compute Stiffness analysis.\n";
		`cp $workdir/$dirCurves/* . &> curvesCopy.log`;
		`cp -r $workdir/$dirCurves/FORCE_CTES .`;
		chmod 0774, "FORCE_CTES";

        	($seq1,$seq2,$firsts,$lasts) = checkPDB($pdb);
	       	print "SEQ1: $seq1\n";
        	print "SEQ2: $seq2\n";
	}
	else{

		print "Curves Analysis has not already been computed, running Curves to obtain necessary files to compute Stiffness analysis.\n";
		mkdir("FORCE_CTES")  if(! -e "FORCE_CTES");
		chmod 0774, "FORCE_CTES";

	        ($seq1,$seq2,$firsts,$lasts) = checkPDB($pdb);
       		print "SEQ1: $seq1\n";
	        print "SEQ2: $seq2\n";

		my $naType = nucleicType($pdb);

	        print "NA_Type: $naType\n";

		###########  CURVES+  #############

		my $iniSnap = 0;
		my $endSnap = 0;
		my $jumpSnap = 1;

	        my $l1 = length($seq1);
        	my $l2 = length($seq2);

	        open NUC,">seq.info";
        	print NUC "  Strand  1 has  $l1 bases (5'-3'): $seq1\n";
	        print NUC "  Strand  2 has  $l2 bases (3'-5'): $seq2\n";
        	print NUC "  NucType: $naType\n";
	        close NUC;

        	my $fullSeq = $seq1.$seq2;

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

		if(-e "$outputCanal.lis"){
			print "WARNING: Saving backup of old lis file: $outputCanal.lis to canalBackup.lis...\n";
			`mv $outputCanal.lis canalBackup.lis`;
		}

		if(-e "$outputCanal.lis"){
		        print "WARNING: Saving backup of old lis file: $outputCanal.lis to canalBackup.lis...\n";
	        	`mv $outputCanal.lis canalBackup.lis`;
		}

		print "Executing Canal program (analysis of Curves+ output)...\n\n";
		#print "Going to execute: exeCanal-histo+series ($outputCanal,$outputCurves,$seq1,1,1,0)\n";
		exeCanal($outputCanal,$outputCurves,$seq1,1,1,0);

		# Saving canal.in/log
		`mv canal.in canalHistoSeries.in`;
		`mv canal.log canalHistoSeries.log`;

		if(! -s "$outputCanal.lis"){
			print "Cannot find Canal output file $outputCanal.lis (or is empty). Exiting...\n";
			exit;
		}
	}

	###########  CANAL Analysis (Nacho Scripts) - stiffness  #############

	# Analisis de stiffness.

	print "Getting Stiffness Constants...\n\n";
	#print "Going to exec force_ctes.pl $seq1 $seq2\n";
	`perl $scriptsDir/force_ctes.pl $seq1 $seq2 &> stiffness.log`;

	#print "Going to exec stiffnessPlots.pl\n";
	chdir("FORCE_CTES");
	`perl $scriptsDir/stiffnessAvgPlots.pl &> stiffnessAvgPlots.log`;
	`perl $scriptsDir/stiffnessTimePlots.pl &> stiffnessTimePlots.log`;

	print "Generating Stiffness HTML tables...\n\n";
	foreach my $file (`ls --color=never *.cte`){
		chomp($file);
		`perl $scriptsDir/stiffnessTableHtml.pl $file $file.html`;
	}
	chdir("../");

	###########  Plotting (Gnuplot Scripts) - Avg Plots  #############

	print "Building Stiffness Constants Plots...\n\n";
	foreach my $file (`ls --color=never FORCE_CTES/*_avg.dat`){
		chomp($file);
		#print "Building Avg plots, file: $file\n";
		if (-s $file){
			my @arr = split '/',$file;
			my $param = $arr[$#arr];
			$param =~s/\_avg\.dat//g;
			$param = ucfirst($param);

			my $xlabel = "Base Pair Step";
			my $units = "Kcal/mol*Degree^2";
			$units = "Kcal/mol*Angstrom^2" if $translation{$param};

			my $ylabel = "$param ($units)";
			my $title = "Helical Parameters Stiffness: $param";

			#print "Command: perl $scriptsDir/plotHelicalParamsStiffness.pl $file $units $stiffnessTable\n";
			#`perl $scriptsDir/plotXYgnuplot.pl $file "$title" "$xlabel" "$ylabel"`;
			#`perl $scriptsDir/plotXYWithSeqgnuplot.pl $file "$title" "$xlabel" "$ylabel"`;
			print "perl $scriptsDir/plotHelicalParamsStiffness.pl $file $units $stiffnessTable &> gnuplot.errors\n";
			`perl $scriptsDir/plotHelicalParamsStiffness.pl $file $units $stiffnessTable &> gnuplot.errors`;
		}
	}

        ###########  Plotting (R Scripts) - Time Plots  #############

	print "Building Time Plots (with associated histograms)...\n\n";
        foreach my $file (`ls --color=never FORCE_CTES/*/*.dat`){
                chomp($file);
		#next if ($file =~ /avg/);
                #print "Building Time plots, file: $file\n";
                if (-s $file){
                        my @arr = split '/',$file;
                        my $param = $arr[$#arr];
                        @arr = split '\.',$param;
                        $param = $arr[0];
                        $param = ucfirst($param);

                        my $xlabel = "Time (Snapshots)";
			my $units = "Degrees";
			$units = "Angstroms" if ($translation{$param});
                        #my $units = "Kcal/mol*Degree^2";
			#$units = "Kcal/mol*Angstrom^2" if $translation{$param};

                        my $ylabel = "$param ($units)";
			my $title = "Helical Parameters  (used in computing Stiffness constants): $param";

                        #print "Command: perl $scriptsDir/plotXYWithSeqgnuplot.pl $file title + x + y \n";
                        #`perl $scriptsDir/plotXYgnuplot.pl $file "$title" "$xlabel" "$ylabel"`;
                        `perl $scriptsDir/plotXY_R.pl $file "$title" "$xlabel" "$ylabel" &> plotR.errors`;
                }
        }

	###########  Building Packed Curves/Canal Serial & Histogram Output files  #############

	my $logTarSer = `tar -zcvf $outputCanal.ser.tgz *.ser`;
	my $logTarHis = `tar -zcvf $outputCanal.his.tgz *.his`;

	###########  Building Output File  #############

	#my $curvesOut = `cat $outputCurves.lis`;
	#my $canalOut = `cat $outputCanal.lis`;
	my $curvesOut = `cat *_curvesOut.lis`;
	my $canalOut = `cat *_canalOut.lis`;
	open OUT, ">nucleicAnalysis.out";
	print OUT "# Curves+ Output File:\n\n";
	print OUT "$curvesOut\n";
	print OUT "# Canal Output File:\n\n";
	print OUT "$canalOut\n";
	close OUT;

	# Building output pdfs (MuG)
	#chdir("FORCE_CTES");
        #my @htmls = `ls --color=never *.html`;
        #my $cathtmls = '';
        #foreach my $html (sort {&byNum} @htmls){
        #        chomp($html);
        #        $cathtmls.= "$html ";
        #}
        #print "wkhtmltopdf.sh --default-header --header-left 'Stiffness Force Constants' --header-spacing 10 toc $cathtmls patata.pdf \n";
        #`wkhtmltopdf.sh --default-header --header-left "Stiffness Force Constants" --header-spacing 10 toc $cathtmls ../../STIFFNESS.ForceConstantsTables.pdf`;

	#`convert *.png ../../STIFFNESS.ForceConstants.pdf`;

	print "DONE!!\n";
}

sub byNum {
	# cgcg.8.cte.html
	my ($tet1,$num1,$tag11,$tag12) = split '\.',$a;
	my ($tet2,$num2,$tag21,$tag22) = split '\.',$b;
	return $num1 <=> $num2;
}
return 1;
