#!/usr/bin/perl -w

use strict;

# Script Input Parameters
my $long=@ARGV;
if ($long<2 or $long >3){
    print "Usage: perl $0 <method> <all-atom> [<system>]\n";
    print "\tWhere method can be: 1 (mcdna), 2 (circular), 3 (prot-dna)\n";
    print "\tWhere all-atom can be: 0 (coarse-grained) or 1 (all-atom)\n";
    print "\tWhere system can be: 0 (structure), 1 (ensemble/trajectory), 3 (both structure and ensemble/trajectory)\n";
    print "Example: perl $0 1 1 0 # To analyse an MCDNA single structure (AA)\n";
    print "Example: perl $0 2 1 1 # To analyse a circular all-atom trajectory\n";
    print "Example: perl $0 3 0 1 # To analyse a Prot-DNA Coarse-Grained trajectory\n";
    print "Example: perl $0 3 1 2 # To analyse a Prot-DNA (AA) structure + trajectory\n";
    exit(0);
}

# Input Parameters.
my ($method,$rebuilt,$traj)=@ARGV;

my %methodTXT;
$methodTXT{"1"} = "mcdna";
$methodTXT{"2"} = "circular";
$methodTXT{"3"} = "prot-dna";

$traj = 0 if !(defined $traj);

my $system = "STRUCT";
$system = "TRAJ" if ($traj == 1);
$system = "STRUCT_AND_TRAJ" if ($traj == 2);

my $resolution = "AA";
$resolution = "CG" if (!$rebuilt);

# Initializing output folders.
my $curr_folder = `pwd`;
chomp($curr_folder);

# Input Sequence 
my $length = 0;
my $seq = '';
my $inputSeq = "$curr_folder/inputSequence.txt";
if (-s "$inputSeq"){
	$length = `wc -c $inputSeq | awk '{print \$1}'`;
	$length += 0;
	$seq = `cat $inputSeq`;
}

# Analysis Scripts
#my $naflex = $ENV{'MCDNA_SCRIPTS'}."/NAFlex";
#my $bending = $ENV{'MCDNA_SCRIPTS'}."/Bending";
#my $elastic = $ENV{'MCDNA_SCRIPTS'}."/ElasticEnergy";
#my $sasa = $ENV{'MCDNA_SCRIPTS'}."/SASA";

my $naflex = "/app/Scripts/MCDNA/NAFlex";
my $bending = "/app/Scripts/MCDNA/Bending";
my $elastic = "/app/Scripts/MCDNA/ElasticEnergy";
my $sasa = "/app/Scripts/MCDNA/SASA";

# Maximum number of atoms for PCAZIP
#my $pcaMaxNatoms = 2000;
my $pcaMaxNatoms = 2500;

# Maximum number of bases for CURVES
my $curvesMaxBases = 250;

print "\n#\n# RESOLUTION: $resolution, SYSTEM: $system, METHOD: $methodTXT{$method}, SEQ_LEN: $length\n#\n\n";

print "STARTING ANALYSIS FOR $methodTXT{$method} METHOD...\n\n"; 

#print "## SECTION NAME ANALYSIS FOR $methodTXT{$method} METHOD\n";
#print "## TRAJ\n";

# MCDNA "regular" METHOD
if ($method == 1){

	if($traj == 1 or $traj == 2){

		print "## TRAJECTORY ANALYSES\n";

		my $out_folder = "$curr_folder/TRAJ_$resolution";

		if (! -s "$out_folder") { print "Ups, $out_folder not found... \nPlease check your inputs and try again!\n\n"; exit;}

		mkdir("$out_folder/ANALYSIS") if (! -s "$out_folder/ANALYSIS");
		chdir("$out_folder/ANALYSIS");

		my $traj_folder = "$out_folder/output_schnarp/display";

		if($rebuilt) {
			# STEP 1: Analysis on Flexibility (NAFlex) 
			print "# STEP 1: Analysis on Flexibility (NAFlex)...\n";
			`perl $naflex/prepTOP.pl $traj_folder/traj.pdb > $traj_folder/struc.prmtop`;
			if ($length < $curvesMaxBases){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Stiffness > NAFlex.stiffness.log 2>&1`;
			}

			`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# End 2 End distances
			if (-s "$out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat"){
				mkdir("$out_folder/ANALYSIS/NAFlex/END-TO-END") if (! -s "$out_folder/ANALYSIS/NAFlex/END-TO-END");
				`cp $out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat $out_folder/ANALYSIS/NAFlex/END-TO-END/distances.dat`;
			}

			my $nats = `grep "^ATOM" $traj_folder/traj.pdb | wc -l`;
			chomp($nats);
			if ($nats < $pcaMaxNatoms){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Pcazip > NAFlex.pcazip.log 2>&1`;
			}

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			#my $numStructs = 10; # HARDCODED!!!!! Need to be find out!!!
			#`Rscript $bending/MuG_DNA_bending_ensemble.R $numStructs $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			print "# STEP 2: Analysis on Bending...\n";
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 0 > elasticEnergy.R.log 2>&1`;
		}
		else{
			# STEP 1: Analysis on Flexibility (NAFlex) 
			print "# STEP 1: Analysis on Flexibility (NAFlex)...\n";
			`perl $naflex/prepTOP.pl $traj_folder/traj.pdb > $traj_folder/struc.prmtop`;
			`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# End 2 End distances
			if (-s "$out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat"){
				mkdir("$out_folder/ANALYSIS/NAFlex/END-TO-END") if (! -s "$out_folder/ANALYSIS/NAFlex/END-TO-END");
				`cp $out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat $out_folder/ANALYSIS/NAFlex/END-TO-END/distances.dat`;
			}

			my $nats = `grep "^ATOM" $traj_folder/traj.pdb | wc -l`;
			chomp($nats);
			if ($nats < $pcaMaxNatoms){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Pcazip > NAFlex.pcazip.log 2>&1`;
			}

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			#my $numStructs = 10; # HARDCODED!!!!! Need to be find out!!!
			#`Rscript $bending/MuG_DNA_bending_ensemble.R $numStructs $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			print "# STEP 2: Analysis on Bending...\n";
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 0 > elasticEnergy.R.log 2>&1`;
		}
	}
	if($traj == 2 or $traj == 0){

		print "## STRUCTURE ANALYSES\n";
 
		my $out_folder = "$curr_folder/EQ_$resolution";

		if (! -s "$out_folder") { print "Ups, $out_folder not found... \nPlease check your inputs and try again!\n\n"; exit;}

		mkdir("$out_folder/ANALYSIS") if (! -s "$out_folder/ANALYSIS");
		chdir("$out_folder/ANALYSIS");

		my $struct_folder = "$out_folder/output_schnarp";

		if ($rebuilt){
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			if ($length < $curvesMaxBases){
				`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
			}

			`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#`Rscript $bending/MuG_DNA_bending_single_structure.R $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 0 > elasticEnergy.R.log 2>&1`;
		}
		else{
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#`Rscript $bending/MuG_DNA_bending_single_structure.R $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 0 > elasticEnergy.R.log 2>&1`;
		}
	}
}
# MCDNA circular METHOD
elsif ($method == 2){  

	if($traj == 1 or $traj == 2){

		print "## TRAJECTORY ANALYSES\n";
		 
		my $out_folder = "$curr_folder/TRAJ_$resolution";

		if (! -s "$out_folder") { print "Ups, $out_folder not found... \nPlease check your inputs and try again!\n\n"; exit;}

		mkdir("$out_folder/ANALYSIS") if (! -s "$out_folder/ANALYSIS");
		chdir("$out_folder/ANALYSIS");

		my $traj_folder = "$out_folder/output_schnarp/display";

		if($rebuilt) {
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			`perl $naflex/prepTOP.pl $traj_folder/traj.pdb > $traj_folder/struc.prmtop`;
			if ($length < $curvesMaxBases){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Stiffness > NAFlex.stiffness.log 2>&1`;
			}

			`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			my $nats = `grep "^ATOM" $traj_folder/traj.pdb | wc -l`;
			chomp($nats);
			if ($nats < $pcaMaxNatoms){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex Pcazip > NAFlex.pcazip.log 2>&1`;
			}

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#my $numStructs = 10; # HARDCODED!!!!! Need to be find out!!!
			#`Rscript $bending/MuG_DNA_bending_ensemble.R $numStructs $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Circularity
			mkdir("$out_folder/ANALYSIS/Circular") if (! -s "$out_folder/ANALYSIS/Circular");
			print "# STEP 3: Analysis on Circularity...\n";
			`Rscript $bending/circle_analysis_webserver_save_csv.R $out_folder $out_folder/ANALYSIS/Circular > Circle.R.log 2>&1`;

			# STEP 4: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 4: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 1 > elasticEnergy.R.log 2>&1`;
		}
		else{
			# STEP 1: Analysis on Flexibility (NAFlex) 
			print "# STEP 1: Analysis on Flexibility (NAFlex)...\n";
			`perl $naflex/prepTOP.pl $traj_folder/traj.pdb > $traj_folder/struc.prmtop`;
			`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			my $nats = `grep "^ATOM" $traj_folder/traj.pdb | wc -l`;
			chomp($nats);
			if ($nats < $pcaMaxNatoms){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex Pcazip > NAFlex.pcazip.log 2>&1`;
			}

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#my $numStructs = 10; # HARDCODED!!!!! Need to be find out!!!
			#`Rscript $bending/MuG_DNA_bending_ensemble.R $numStructs $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Circularity
			mkdir("$out_folder/ANALYSIS/Circular") if (! -s "$out_folder/ANALYSIS/Circular");
			print "# STEP 3: Analysis on Circularity...\n";
			`Rscript $bending/circle_analysis_webserver_save_csv.R $out_folder $out_folder/ANALYSIS/Circular > Circle.R.log 2>&1`;

			# STEP 4: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 4: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 1 > elasticEnergy.R.log 2>&1`;
		}
	}
	if($traj == 2 or $traj == 0){

		print "## STRUCTURE ANALYSES\n";

		my $out_folder = "$curr_folder/EQ_$resolution";

		if (! -s "$out_folder") { print "Ups, $out_folder not found... \nPlease check your inputs and try again!\n\n"; exit;}

		mkdir("$out_folder/ANALYSIS") if (! -s "$out_folder/ANALYSIS");
		chdir("$out_folder/ANALYSIS");

		my $struct_folder = "$out_folder/output_schnarp";

		if ($rebuilt){
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			if ($length < $curvesMaxBases){
				`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
			}
			`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#`Rscript $bending/MuG_DNA_bending_single_structure.R $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Circularity
			mkdir("$out_folder/ANALYSIS/Circular") if (! -s "$out_folder/ANALYSIS/Circular");
			print "# STEP 3: Analysis on Circularity...\n";
			`Rscript $bending/circle_analysis_webserver_save_csv.R $out_folder $out_folder/ANALYSIS/Circular > Circle.R.log 2>&1`;

			# STEP 4: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 4: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 1 > elasticEnergy.R.log 2>&1`;
		}
		else{
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#`Rscript $bending/MuG_DNA_bending_single_structure.R $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Circularity
			mkdir("$out_folder/ANALYSIS/Circular") if (! -s "$out_folder/ANALYSIS/Circular");
			print "# STEP 3: Analysis on Circularity...\n";
			`Rscript $bending/circle_analysis_webserver_save_csv.R $out_folder $out_folder/ANALYSIS/Circular > Circle.R.log 2>&1`;

			# STEP 4: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 4: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 0 1 > elasticEnergy.R.log 2>&1`;
		}
	}
}
# MCDNA Prot-DNA METHOD
elsif ($method == 3){

	if($traj == 1 or $traj == 2){

		print "## TRAJECTORY ANALYSES\n";

		my $out_folder = "$curr_folder/TRAJ_$resolution";

		if (! -s "$out_folder") { print "Ups, $out_folder not found... \nPlease check your inputs and try again!\n\n"; exit;}

		mkdir("$out_folder/ANALYSIS") if (! -s "$out_folder/ANALYSIS");
		chdir("$out_folder/ANALYSIS");

		my $traj_folder = "$out_folder/output_schnarp/display";

		if($rebuilt) {
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			#`perl $naflex/prepTOP.pl $traj_folder/traj.pdb > $traj_folder/struc.prmtop`;

			#`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
			if ($length < $curvesMaxBases){
				#`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/nucleic_mc_dna_str.pdb $traj_folder/struc.prmtop $traj_folder/nucleic_mc_dna_str.dcd $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/nucleic_mc_dna_str.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.nc $out_folder/ANALYSIS NAFlex Stiffness > NAFlex.stiffness.log 2>&1`;
			}

			`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# End 2 End distances
			if (-s "$out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat"){
				mkdir("$out_folder/ANALYSIS/NAFlex/END-TO-END") if (! -s "$out_folder/ANALYSIS/NAFlex/END-TO-END");
				`cp $out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat $out_folder/ANALYSIS/NAFlex/END-TO-END/distances.dat`;
			}

			my $nats = `grep "^ATOM" $traj_folder/traj.pdb | wc -l`;
			chomp($nats);
			if ($nats < $pcaMaxNatoms){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex Pcazip > NAFlex.pcazip.log 2>&1`;
			}

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#my $numStructs = 10; # HARDCODED!!!!! Need to be find out!!!
			#`Rscript $bending/MuG_DNA_bending_ensemble.R $numStructs $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 1 0 > elasticEnergy.R.log 2>&1`;

			# STEP 4: Analysis on SASA (Virtual Footprinting)
			mkdir("$out_folder/ANALYSIS/Sasa") if (! -s "$out_folder/ANALYSIS/Sasa");
			print "# STEP 4: Analysis on SASA (Virtual Footprinting)...\n";
			chdir("$out_folder/ANALYSIS/Sasa");
			`perl $naflex/prepTOP.pl $traj_folder/traj.pdb > $traj_folder/struc.proteins.prmtop`;
			#`perl $sasa/getSASA.pl $traj_folder/traj.pdb $traj_folder/struc.proteins.prmtop $traj_folder/traj.dcd > SASA.log 2>&1`;
			`perl $sasa/getSASA.pl $traj_folder/traj.pdb $traj_folder/traj.pdb $traj_folder/traj.dcd > SASA.log 2>&1`;
			chdir("..");
			
		}
		else{
			# STEP 1: Analysis on Flexibility (NAFlex) 
			print "# STEP 1: Analysis on Flexibility (NAFlex)...\n";
			`perl $naflex/prepTOP.pl $traj_folder/traj.pdb > $traj_folder/struc.prmtop`;
			`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# End 2 End distances
			if (-s "$out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat"){
				mkdir("$out_folder/ANALYSIS/NAFlex/END-TO-END") if (! -s "$out_folder/ANALYSIS/NAFlex/END-TO-END");
				`cp $out_folder/ANALYSIS/NAFlex/CONTACTS/NUC-NUC/end2end.dat $out_folder/ANALYSIS/NAFlex/END-TO-END/distances.dat`;
			}

			my $nats = `grep "^ATOM" $traj_folder/traj.pdb | wc -l`;
			chomp($nats);
			if ($nats < $pcaMaxNatoms){
				`perl $naflex/runNucleicAcidAnalysis.pl $traj_folder/traj.pdb $traj_folder/struc.prmtop $traj_folder/traj.dcd $out_folder/ANALYSIS NAFlex Pcazip > NAFlex.pcazip.log 2>&1`;
			}

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#my $numStructs = 10; # HARDCODED!!!!! Need to be find out!!!
			#`Rscript $bending/MuG_DNA_bending_ensemble.R $numStructs $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 1 0 > elasticEnergy.R.log 2>&1`;

			# STEP 4: Analysis on SASA (Virtual Footprinting)
			mkdir("$out_folder/ANALYSIS/Sasa") if (! -s "$out_folder/ANALYSIS/Sasa");
			print "# STEP 4: Analysis on SASA (Virtual Footprinting)...\n";
			chdir("$out_folder/ANALYSIS/Sasa");
			`perl $sasa/getSASA.pl $traj_folder/traj.pdb $traj_folder/traj.pdb $traj_folder/traj.dcd > SASA.log 2>&1`;
			chdir("..");
		}
	}
	if($traj == 2 or $traj == 0){

		print "## STRUCTURE ANALYSES\n";

		my $out_folder = "$curr_folder/EQ_$resolution";

		if (! -s "$out_folder") { print "Ups, $out_folder not found... \nPlease check your inputs and try again!\n\n"; exit;}

		mkdir("$out_folder/ANALYSIS") if (! -s "$out_folder/ANALYSIS");
		chdir("$out_folder/ANALYSIS");

		my $struct_folder = "$out_folder/output_schnarp";

		if ($rebuilt){
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			if ($length < $curvesMaxBases){
				`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.curves.pdb $out_folder/ANALYSIS NAFlex Curves > NAFlex.curves.log 2>&1`;
			}
			`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#`Rscript $bending/MuG_DNA_bending_single_structure.R $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 1 0 > elasticEnergy.R.log 2>&1`;

			# STEP 4: Analysis on SASA (Virtual Footprinting)
			mkdir("$out_folder/ANALYSIS/Sasa") if (! -s "$out_folder/ANALYSIS/Sasa");
			print "# STEP 4: Analysis on SASA (Virtual Footprinting)...\n";
			chdir("$out_folder/ANALYSIS/Sasa");
			`perl $naflex/prepTOP.pl $struct_folder/str.pdb > $struct_folder/str.top`;
			`perl $sasa/getSASA_Pdb.pl $struct_folder/str.pdb > SASA.log 2>&1`;
			chdir("..");
			
		}
		else{
			# STEP 1: Analysis on Flexibility (NAFlex)... 
			print "# STEP 1: Analysis on Flexibility (NAFlex)\n";
			`perl $naflex/runNucleicAcidAnalysisPDB.pl $struct_folder/str.pdb $out_folder/ANALYSIS NAFlex DistanceContactMaps > NAFlex.dist.log 2>&1`;

			# STEP 2: Analysis on Bending
			mkdir("$out_folder/ANALYSIS/Bending") if (! -s "$out_folder/ANALYSIS/Bending");
			print "# STEP 2: Analysis on Bending...\n";
			#`Rscript $bending/MuG_DNA_bending_single_structure.R $out_folder $out_folder/ANALYSIS/Bending > Bending.R.log 2>&1`;
			`Rscript $bending/MuG_DNA_bending_extended_save_csv.R $out_folder $out_folder/ANALYSIS/Bending > newBending.R.log 2>&1`;

			# STEP 3: Analysis on Elastic Energy
			mkdir("$out_folder/ANALYSIS/ElasticEnergy") if (! -s "$out_folder/ANALYSIS/ElasticEnergy");
			print "# STEP 3: Analysis on Elastic Energy...\n";
			`Rscript $elastic/MCDNA_comp_elastic_ene_mcdna.R $inputSeq $out_folder $out_folder/ANALYSIS/ElasticEnergy 1 0 > elasticEnergy.R.log 2>&1`;

			# STEP 4: Analysis on SASA (Virtual Footprinting)
			mkdir("$out_folder/ANALYSIS/Sasa") if (! -s "$out_folder/ANALYSIS/Sasa");
			print "# STEP 4: Analysis on SASA (Virtual Footprinting)...\n";
			chdir("$out_folder/ANALYSIS/Sasa");
			`perl $sasa/getSASA_Pdb.pl $struct_folder/str.pdb > SASA.log 2>&1`;
			chdir("..");
			
		}
	}
}

# Building downloadable tar package
print "Building downloadable tar package for MCDNA Analysis...\n";
chdir("$curr_folder");
#`tar -zcvf mcdna.analysis.tgz EQ*/ANALYSIS` if ($traj == 0);
#`tar -zcvf mcdna.analysis.tgz TRAJ*/ANALYSIS` if ($traj == 1);
#`tar -zcvf mcdna.analysis.tgz EQ*/ANALYSIS TRAJ*/ANALYSIS` if ($traj == 2);

&download_readme($seq,$methodTXT{$method},$resolution,$system);
&download_eq($methodTXT{$method},$resolution);
&download_traj($methodTXT{$method},$resolution);

chdir("$curr_folder/download");
my $m = $methodTXT{$method};
my $tf = "$m-$resolution-$system".".tgz";
`tar -zcvf $tf *`;
chdir("..");
`mv download/$tf .`;

# Cleaning downloadable (and duplicated) info
#`rm -r download`;

print "DONE!!\n";

################################################

sub download_readme {
	my ($seq,$method,$resolution,$system) = @_;

	mkdir("download") if (! -s "download");

	my $nstructs = 0;
	my $log = "TRAJ_$resolution/output_schnarp/display/cpptraj_dcd.log";
	if (-s "$log"){
		my $g = `grep "frames and processed" $log`;
		$g =~ /Read (\d+) frames and processed/;
		$nstructs = $1;
	}

	# Auxiliary Files
	open SUM,">download/summary.txt";
	print SUM "# MCDNA process summary #\n";
	print SUM "Sequence: $seq\n";
	print SUM "Method: $method\n";
	print SUM "Resolution: $resolution\n";
	print SUM "System: $system\n";
	print SUM "Number of structures: $nstructs\n" if ($system =~ /TRAJ/ and $nstructs);
	close SUM;	
}

sub download_eq {
	my ($method,$resolution) = @_;

	mkdir("download") if (! -s "download");
	
	# MCDNA Analysis
	mkdir("download/STRUCT") if (! -s "download/STRUCT");
	mkdir("download/STRUCT/ANALYSIS") if (! -s "download/STRUCT/ANALYSIS");

	`cp -r EQ_$resolution/ANALYSIS/Bending download/STRUCT/ANALYSIS`;
	`cp -r EQ_$resolution/ANALYSIS/ElasticEnergy download/STRUCT/ANALYSIS`;

	mkdir("download/STRUCT/ANALYSIS/Contacts") if (! -s "download/STRUCT/ANALYSIS/Contacts");
	`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CONTACTS/*/*.dat download/STRUCT/ANALYSIS/Contacts`;
	`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CONTACTS/*/*.png download/STRUCT/ANALYSIS/Contacts`;

	if ($resolution eq "AA"){
		mkdir("download/STRUCT/ANALYSIS/HelicalParams") if (! -s "download/STRUCT/ANALYSIS/HelicalParams");
		mkdir("download/STRUCT/ANALYSIS/HelicalParams/grooves") if (! -s "download/STRUCT/ANALYSIS/HelicalParams/grooves");
		mkdir("download/STRUCT/ANALYSIS/HelicalParams/helical_bp") if (! -s "download/STRUCT/ANALYSIS/HelicalParams/helical_bp");
		mkdir("download/STRUCT/ANALYSIS/HelicalParams/helical_bps") if (! -s "download/STRUCT/ANALYSIS/HelicalParams/helical_bps");
		mkdir("download/STRUCT/ANALYSIS/HelicalParams/axis_bp") if (! -s "download/STRUCT/ANALYSIS/HelicalParams/axis_bp");
		mkdir("download/STRUCT/ANALYSIS/HelicalParams/backbone_torsions") if (! -s "download/STRUCT/ANALYSIS/HelicalParams/backbone_torsions");
		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/helical_bp/*.png download/STRUCT/ANALYSIS/HelicalParams/helical_bp`;
		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/helical_bp/*.dat download/STRUCT/ANALYSIS/HelicalParams/helical_bp`;

		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/helical_bpstep/*.png download/STRUCT/ANALYSIS/HelicalParams/helical_bps`;
		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/helical_bpstep/*.dat download/STRUCT/ANALYSIS/HelicalParams/helical_bps`;

		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/grooves/*.png download/STRUCT/ANALYSIS/HelicalParams/grooves`;
		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/grooves/*.dat download/STRUCT/ANALYSIS/HelicalParams/grooves`;

		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/axis_bp/*.png download/STRUCT/ANALYSIS/HelicalParams/axis_bp`;
		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/axis_bp/*.dat download/STRUCT/ANALYSIS/HelicalParams/axis_bp`;

		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/backbone_torsions/*.png download/STRUCT/ANALYSIS/HelicalParams/backbone_torsions`;
		`cp EQ_$resolution/ANALYSIS/NAFlex/PDB/CURVES/backbone_torsions/*.dat download/STRUCT/ANALYSIS/HelicalParams/backbone_torsions`;

	}
	if ($method eq "circular"){
		`cp -r EQ_$resolution/ANALYSIS/Circular download/STRUCT/ANALYSIS`;
	}

	if ($method eq "prot-dna"){
		mkdir("download/STRUCT/ANALYSIS/Sasa") if (! -s "download/STRUCT/ANALYSIS/Sasa");
		`cp EQ_$resolution/ANALYSIS/Sasa/sasa.*.dat download/STRUCT/ANALYSIS/Sasa`;
	}
}

sub download_traj {
	my ($method,$resolution) = @_;

	mkdir("download") if (! -s "download");
	
	# MCDNA Analysis
	mkdir("download/TRAJ") if (! -s "download/TRAJ");
	mkdir("download/TRAJ/ANALYSIS") if (! -s "download/TRAJ/ANALYSIS");

	`cp -r TRAJ_$resolution/ANALYSIS/Bending download/TRAJ/ANALYSIS`;
	`cp -r TRAJ_$resolution/ANALYSIS/ElasticEnergy download/TRAJ/ANALYSIS`;

	mkdir("download/TRAJ/ANALYSIS/Contacts") if (! -s "download/TRAJ/ANALYSIS/Contacts");
	`cp TRAJ_$resolution/ANALYSIS/NAFlex/CONTACTS/*/*.dat download/TRAJ/ANALYSIS/Contacts`;
	`cp TRAJ_$resolution/ANALYSIS/NAFlex/CONTACTS/*/*.png download/TRAJ/ANALYSIS/Contacts`;

	if (-s "TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP"){
		mkdir("download/TRAJ/ANALYSIS/Pca") if (! -s "download/TRAJ/ANALYSIS/Pca");
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP/*.pdb download/TRAJ/ANALYSIS/Pca`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP/*.dat download/TRAJ/ANALYSIS/Pca`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP/*.dat.png download/TRAJ/ANALYSIS/Pca`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP/NAFlex_pcazipOut.bfactors download/TRAJ/ANALYSIS/Pca`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP/NAFlex_pcazipOut.collectivity download/TRAJ/ANALYSIS/Pca`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP/NAFlex_pcazipOut.evals download/TRAJ/ANALYSIS/Pca`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/PCAZIP/NAFlex_pcazipOut.stiffness download/TRAJ/ANALYSIS/Pca`;
	}

	mkdir("download/TRAJ/ANALYSIS/EndToEnd") if (! -s "download/TRAJ/ANALYSIS/EndToEnd");
	`cp TRAJ_$resolution/ANALYSIS/NAFlex/END-TO-END/distances.dat download/TRAJ/ANALYSIS/EndToEnd`;

	#mkdir("download/TRAJ/ANALYSIS/Stiffness") if (! -s "download/TRAJ/ANALYSIS/Stiffness");
	#`cp TRAJ_$resolution/ANALYSIS/NAFlex/STIFFNESS/* download/TRAJ/ANALYSIS/Stiffness`;

	if ($resolution eq "AA"){
		mkdir("download/TRAJ/ANALYSIS/HelicalParams") if (! -s "download/TRAJ/ANALYSIS/HelicalParams");
		mkdir("download/TRAJ/ANALYSIS/HelicalParams/grooves") if (! -s "download/TRAJ/ANALYSIS/HelicalParams/grooves");
		mkdir("download/TRAJ/ANALYSIS/HelicalParams/helical_bp") if (! -s "download/TRAJ/ANALYSIS/HelicalParams/helical_bp");
		mkdir("download/TRAJ/ANALYSIS/HelicalParams/helical_bps") if (! -s "download/TRAJ/ANALYSIS/HelicalParams/helical_bps");
		mkdir("download/TRAJ/ANALYSIS/HelicalParams/axis_bp") if (! -s "download/TRAJ/ANALYSIS/HelicalParams/axis_bp");
		mkdir("download/TRAJ/ANALYSIS/HelicalParams/backbone_torsions") if (! -s "download/TRAJ/ANALYSIS/HelicalParams/backbone_torsions");
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/helical_bp/*.png download/TRAJ/ANALYSIS/HelicalParams/helical_bp`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/helical_bp/*.dat download/TRAJ/ANALYSIS/HelicalParams/helical_bp`;

		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/helical_bpstep/*.png download/TRAJ/ANALYSIS/HelicalParams/helical_bps`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/helical_bpstep/*.dat download/TRAJ/ANALYSIS/HelicalParams/helical_bps`;

		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/grooves/*.png download/TRAJ/ANALYSIS/HelicalParams/grooves`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/grooves/*.dat download/TRAJ/ANALYSIS/HelicalParams/grooves`;

		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/axis_bp/*.png download/TRAJ/ANALYSIS/HelicalParams/axis_bp`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/axis_bp/*.dat download/TRAJ/ANALYSIS/HelicalParams/axis_bp`;

		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/backbone_torsions/*.png download/TRAJ/ANALYSIS/HelicalParams/backbone_torsions`;
		`cp TRAJ_$resolution/ANALYSIS/NAFlex/CURVES/backbone_torsions/*.dat download/TRAJ/ANALYSIS/HelicalParams/backbone_torsions`;

	}
	if ($method eq "circular"){
		`cp -r EQ_$resolution/ANALYSIS/Circular download/STRUCT/ANALYSIS`;
	}
	if ($method eq "prot-dna"){
		mkdir("download/TRAJ/ANALYSIS/Sasa") if (! -s "download/TRAJ/ANALYSIS/Sasa");
		`cp TRAJ_$resolution/ANALYSIS/Sasa/sasa.*.dat download/TRAJ/ANALYSIS/Sasa`;
	}
}

