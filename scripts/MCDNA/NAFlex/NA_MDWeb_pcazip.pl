#!/usr/bin/perl -w

use strict;

# run_pcazip
# Executing pcazip
# Saving output files in $workdir/Pcazip
sub run_pcazip{

	my ($workdir,$prefix,$pdb,$crd) = @_;

	chdir($workdir);
	#mkdir("Pcazip_$prefix")  if(! -e "Pcazip_$prefix");
	#chdir("Pcazip_$prefix");
	mkdir("PCAZIP")  if(! -e "PCAZIP");
	chdir("PCAZIP");

	my $output = $prefix."_pcazipOut";

	my $nvecs = 10;
	
	# Compressing Trajectory to PCZ file.
	print "Executing PCAzip Analysis, first compressing the trajectory...\n\n";
	#print "Going to execute: exePcazip($pdb,$output,$nvecs)\n";
	exePcazip($pdb,$crd,$output,$nvecs);

	if(! -s "$output.pcz"){
		print "Cannot find Pcazip output file $output.pcz (or is empty). Exiting...\n";
		exit;
	}

	# Computing Bfactors from PCZ.
	print "Computing Bfactors from PCZ...\n\n";
	computeBfactors($output);

	# Computing Trajectory Projections from PCZ.
	print "Computing Trajectory Projections from PCZ...\n";
	foreach (my $i=1; $i<=$nvecs; $i++){
		print "\tProjections to Vector $i...\n";
		computeProjections($output,$i);

		my $title = "Trajectory Projection to Vector $i";
		my $xlabel = "Time (Snapshots)";
		my $ylabel = "Displacement (Angstroms)";
		#print "\tPlotting plotXY_R.pl $output.proj$i.dat $title $xlabel $ylabel &> gnuplot.errors\n";
		#`perl $scriptsDir/plotXY_R.pl $output.proj$i.plot.dat "$title" "$xlabel" "$ylabel" &> gnuplot.errors`;
		`python $scriptsDir/plotXY_Python.py $output.proj$i.plot.dat "$title" "$xlabel" "$ylabel" &> gnuplot.errors`;
	}

	# Computing Animations from PCZ.
	print "Computing Animations from PCZ...\n";
	foreach (my $i=1; $i<=$nvecs; $i++){
		print "\tAnimation for Vector $i...\n";
		computeAnimations($output,$i);
	}

	# Computing PCZ info.
	print "Computing PCZ info: evals & collectivity indexes.\n";
	computePczInfo($output);

	###########  Building Output File  #############

	my $pcaOut = `cat pcazip.log`;
	open OUT, ">nucleicAnalysis.out";
	print OUT "# Pcazip Output File:\n\n";
	print OUT "$pcaOut\n";
	close OUT;

        ###########  MuG output  #############

	#for (my $numAnim=1; $numAnim<=10; $numAnim++){
	#	`cp $output.anim$numAnim.pdb ../PCAZIP.animEvec$numAnim.pdb`;
	#} 

	print "DONE!!\n";
}

###########################  SUBROUTINES  ################################

# exePcazip
# Executing PcatoolS program
sub exePcazip {
	my ($pdb,$crd,$out,$nvecs) = @_;

	# Executing Pcazip
	#`$pcazip -i $crd -p $pdb -q 90 -e $nvecs -o $out &> pcazip.log`;
	print "$pcazip -i $crd -o $out.pcz -q 90 -p $pdb -v -e $nvecs -M ~\@H* \n";
	my $outPca = `$pcazip -i $crd -o $out.pcz -q 90 -p $pdb -v -e $nvecs -M ~\@H* &> pcazip.log`;
	#my $outPcaGauss = `$pcazip -i $crd -o $out.gauss.pcz -q 90 -p $pdb -g -v -e $nvecs -M ~\@H* &> pcazip.gauss.log`;
}

# computeBfactors
# Computing Bfactors with pczdump program from a PCZ file.
sub computeBfactors {
	my ($out) = @_;

	# Executing pczdump
	my $outPczdump = `$pczdump -i $out.pcz --fluc=0 --bfactor -o $out.bfactors &> pcazdump.bfactors.log`;
}

# computePczInfo
# Computing PCZ info with pczdump program from a PCZ file.
sub computePczInfo {
	my ($out) = @_;

	# Executing pczdump
	my $outPczdump;
	$outPczdump = `$pczdump -i $out.pcz --eval -o $out.evals &> pcazdump.evals.log`;
	$outPczdump = `$pczdump -i $out.pcz --forcecte 300 -o $out.stiffness &> pcazdump.stiffness.log`;
	$outPczdump = `$pczdump -i $out.pcz --collectivity -o $out.collectivity &> pcazdump.collectivity.log`;
	$outPczdump = `$pczdump -i $out.pcz --info -o $out.info &> pcazdump.info.log`;
}

# computeProjections
# Computing Trajectory Projections with pczdump program from a PCZ file.
sub computeProjections {
	my ($out,$numProj) = @_;

	# Executing pczdump
	my $outPczdump = `$pczdump -i $out.pcz --proj $numProj > $out.proj$numProj.dat`;

	open PLOT, ">$out.proj$numProj.plot.dat";
	open DAT, "$out.proj$numProj.dat";
	my $i = 1;
	while(<DAT>){
                chomp;
		# -51.374 -54.501-107.734 -6.296 -37.502 -25.682  -4.362 -36.470  42.396 -73.121
		#my @arr = split ' ';
		my $line = $_;
		# Here M is added every 8 characters as a code to easily obtain all projections (split the line in segments -floats- of length 8)
		$line =~ s/(.{8})/$1M/g; 
		# Now M is used to split the line.
		my @arr = split 'M',$line;
	
		foreach (my $j=0;$j<=$#arr;$j++){	
			my $val = $arr[$j];
	                $val =~s/\s+//g;
	                print PLOT "$i $val\n";
        	        $i++;
		}
	}
	close DAT;
	close PLOT;

	my ($mean,$stdev) = getAvg("$out.proj$numProj.plot.dat");
	open STDEV,">$out.proj$numProj.stats";
	print STDEV "Mean: $mean, StDev: $stdev\n";
	close STDEV;
}

# computeAnimations
# Computing Animations with pczdump program from a PCZ file.
sub computeAnimations {
	my ($out,$numProj) = @_;

	# Executing pczdump
	my $outPczdump = `$pczdump -i $out.pcz --anim=$numProj --pdb -o $out.anim$numProj.pdb`;
}

return 1;
