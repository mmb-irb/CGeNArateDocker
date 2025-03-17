#!/usr/bin/perl -w

use strict;

# Script Input Parameters
my $long=@ARGV;
if ($long<5 or $long >6){
    print "Usage: perl $0 <seq_file> <num_structures> <mode> <mcdna_prot config file> <all-atom> [0 (or None): eqStruct, 1: <create_traj>, 2: <create_both_eqStruct_and_traj>]\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 3 mcdna_prot.config 1 0 # To generate a single structure (AA) (equilibrated)\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 3 mcdna_prot.config 1 1 # To generate an all-atom trajectory\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 3 mcdna_prot.config 0 1 # To generate a Coarse-Grained trajectory\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 3 mcdna_prot.config 1 2 # To generate both: (AA) structure + trajectory\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 3 mcdna_prot.config 0 2 # To generate both: (CG) structure + trajectory\n";
    exit(0);
}

# Input Parameters.
my ($seq,$nstructs,$mode,$config,$rebuilt,$traj)=@ARGV;

$traj = 0 if !(defined $traj);

my $system = "STRUCT";
$system = "TRAJ" if ($traj == 1);
$system = "STRUCT_AND_TRAJ" if ($traj == 2);

my $resolution = "AA";
$resolution = "CG" if (!$rebuilt);

# Auxiliar scripts/binaries
my $appFolder = "/app/Scripts/cgenarate-materials/Proteins";
my $scriptsFolder = "/app/Scripts/MCDNA/NAFlex";
my $MCDNAbinStruct_AA = "getCG_AmberAA.sh";

# MCDNA 2025 extended PATH.
$ENV{'PATH'}="$ENV{'PATH'}:/app/Scripts/cgenarate-materials/fdhelix/";

# Input Sequence. Check existence and number of nucleotides (*2 for 2 strands).
if (! -s "$seq"){
	print "Ups, input file $seq not found! Please check it and try again!\n";
	exit;
}

# Sequence length
open FILE,"$seq";
my $arr = <FILE>;
chomp($arr);
close FILE;

# Sequence length 
my $seq_len = length($arr);

# Sequence length *2 (2 strands)
my $num = length($arr);
$num *= 2;

# Number of protein to consider
my $num_prot = `wc -l $config  | cut -f1 -d" "`;
chomp($num_prot);

# Initializing output folders.
my $curr_folder = `pwd`;
chomp($curr_folder);

print "## BUILDING CONFIGURATION FILES\n";

# STEP 1: Building CGeNArate-prot Configuration File
print "# STEP 1: Building CGeNArate-prot Configuration File\n";

my %prots;
my @alphabet = ('C'..'Z');
open PROTINFO,">proteins.info";
my $i = 0;

my $aux_folder = "pdb_inputs";
mkdir("$aux_folder") if (! -s "$aux_folder");

my $configFile = "mcdna_prot.config";
open CONFIG,">$configFile";
open INCONFIG,"$config";
while(<INCONFIG>){
        chomp;
        # 12 1iv6 10
        # 13 1j46 30
        my ($l,$pdb,$ini_pos) = split ' ';
        #`wget http://mmb.irbbarcelona.org/api/pdb/$pdb.pdb`;
        `wget https://www.ebi.ac.uk/pdbe/entry-files/download/pdb$pdb.ent`;

        my $gr = `grep ^MODEL pdb$pdb.ent`;
        if ($gr){
            open OPDB,">$aux_folder/$pdb.pdb";
            open IPDB,"pdb$pdb.ent";
            while(<IPDB>){
		if ($_ =~/ENDMDL/){
		    last;
		}
		print OPDB "$_";
            }
            close IPDB;
            close OPDB;
            `rm pdb$pdb.ent`;
        }
        else{
            `mv pdb$pdb.ent $aux_folder/$pdb.pdb`;
        }
        #print CONFIG "$ini_pos $protdnapath/$pdb.helprms $l\n";
        print CONFIG "$ini_pos $pdb $l\n";

	my $chain = $alphabet[$i];

	# Protein PDB-chain match for the Web NGL visualizer (Contact plots)
	print PROTINFO "$pdb:$chain\n";

	$i++;
}
close INCONFIG;
close CONFIG;
close PROTINFO;

print "# STEP 2: Building structure\n";
`python $appFolder/buildProteins.py`;
	
my $out_folder = "$curr_folder/TRAJ_$resolution";
mkdir("$out_folder") if (! -s "$out_folder");
my $out_folder2 = "$curr_folder/TRAJ_$resolution/output_schnarp";
mkdir("$out_folder2") if (! -s "$out_folder2");
my $out_folder3 = "$curr_folder/TRAJ_$resolution/output_schnarp/display";
mkdir("$out_folder3") if (! -s "$out_folder3");

# Fixing PDB file with chains
`python $appFolder/convertPDB.py input/Proteins.pdb input/ProteinsFixed.pdb`;
`cp input/ProteinsFixed.pdb $out_folder3/traj.pdb`;

print "# STEP 3: Executing CGenerate...\n";
`cp $appFolder/input/*.txt input`;
&runCgenerate($appFolder,$out_folder,$seq,$nstructs,$rebuilt);

# Building downloadable tar package
print "Building downloadable tar package...\n";
chdir("$curr_folder");
#`tar -zcvf mcdna.tgz EQ*` if ($traj == 0);
#`tar -zcvf mcdna.tgz TRAJ*` if ($traj == 1);
#`tar -zcvf mcdna.tgz EQ* TRAJ*` if ($traj == 2);

my $method = "prot-dna";
#&download_readme($arr,$method,$resolution,$system);
#&download_eq($resolution);
&download_traj($resolution);

chdir("$curr_folder/download");
my $tf = "$method-$resolution-$system".".tgz";
`tar -zcvf $tf *`;
chdir("..");
`mv download/$tf .`;

# Cleaning downloadable (and duplicated) info
#`rm -r download`;

print "DONE!!\n";

##################################################################################################################

sub download_traj {
        my ($resolution) = @_;

        mkdir("download") if (! -s "download");

        # CGeNArate Structure
        mkdir("download/TRAJ") if (! -s "download/TRAJ");
        mkdir("download/TRAJ/CGeNArate") if (! -s "download/TRAJ/CGeNArate");

        `cp TRAJ_$resolution/output_schnarp/display/traj.pdb download/TRAJ/CGeNArate`;
        `cp TRAJ_$resolution/output_schnarp/display/traj.dcd download/TRAJ/CGeNArate`;
        #`cp -r TRAJ_$resolution/output_tables_helpar download/TRAJ/CGeNArate`;

        #`cp TRAJ_$resolution/output_schnarp/display/mc_dna_str.pdb download/TRAJ/CGeNArate`;
        #`cp TRAJ_$resolution/output_schnarp/display/mc_dna_str.dcd download/TRAJ/CGeNArate`;
        `cp proteins.info download/TRAJ/CGeNArate`;
        `cp input/*.pdb download/TRAJ/CGeNArate`;
}

sub runCgenerate {

        my ($dir_binary,$out_folder,$seq,$nstructs,$atom) = @_;

        mkdir("$out_folder") if (!-s "$out_folder");
        mkdir("$out_folder/output_pdb") if (!-s "$out_folder/output_pdb");
        mkdir("$out_folder/output_schnarp") if (!-s "$out_folder/output_schnarp");

        # CGeNArate only works if executed where the binary is (uses hard-coded auxiliary files)
        #print "Changing folder to: $dir_binary\n";
        #chdir("$dir_binary");

        my $config = "$out_folder/output_schnarp/cgenerate_config.toml";
        open CONFIG,">$config";
        print CONFIG "# CGeNArate config file, created by CGeNArate web server #\n";
        print CONFIG "[simulation]\n";
        print CONFIG "name = 'CGeNArate web server'\n";
        print CONFIG "T0 = 298.0 # Temperature (K)\n";
        print CONFIG "dt0 = 0.1 # timestep in picoseconds\n";
        print CONFIG "N = $num # Number of beads\n";
        print CONFIG "frames = $nstructs # Number of frames\n";
        print CONFIG "steps = 200 # Number of steps between frames\n";
        print CONFIG "linear = true\n";
        print CONFIG "seed = 1234 # random seed\n\n";
        print CONFIG "[input]\n";
        print CONFIG "file = 'input/Proteins.pdb'\n\n";
        print CONFIG "[output]\n";
        print CONFIG "coordinates = '$out_folder/output_schnarp/cgenarate_traj.mdcrd'\n";
        print CONFIG "energy = '$out_folder/output_schnarp/cgenarate_energy.csv'\n";
        print CONFIG "plots = '$out_folder/output_schnarp/cgenarate_plots.csv'\n\n";
        print CONFIG "[restore]\n";
        print CONFIG "mode = 'none'\n";
        #print CONFIG "current = '$out_folder/output_schnarp/current.pdb'\n";
        #print CONFIG "restart = '$out_folder/output_schnarp/restart.pdb'\n";
        close CONFIG;

        #print "$dir_binary/CGeNArate.seq.exe $config\n";
        #`$dir_binary/CGeNArate.seq.exe $config`;

	print "cat data.in | $appFolder/CGeNArateProteins $config\n";
	`cat data.in | $appFolder/CGeNArateProteins $config`;

        if($atom){
                # Generating GLIMPS reference structure

                #mkdir("$out_folder") if (!-s "$out_folder");
                #mkdir("$out_folder/output_schnarp") if (!-s "$out_folder/output_schnarp");

                my $curr_folder = `pwd`;
                chomp($curr_folder);

                mkdir("$out_folder/output_pdb") if (!-s "$out_folder/output_pdb");
                chdir("$out_folder/output_pdb");
                print("Changing to: $out_folder/output_pdb\n");

		print "$MCDNAbinStruct_AA $seq $out_folder/output_pdb/structure_000000_ \n";
		`$MCDNAbinStruct_AA $seq ${out_folder}/output_pdb/structure_000000_`;
		#`cp $out_folder/output_pdb/structure_000000_CG.pdb $out_folder/output_pdb/structure_000000.pdb`;
		`cp $out_folder/output_pdb/structure_000000_AA_fdh.pdb $out_folder/output_schnarp/reference.pdb`;

                chdir("$curr_folder");

                # Building AA trajectory with GLIMPS
                print "python $appFolder/Rebuild_nmerProteins.py --ref $out_folder/output_schnarp/reference.pdb --inputtop input/ProteinsChains.pdb --inputtraj $out_folder/output_schnarp/cgenarate_traj.mdcrd --proteintop $out_folder/output_schnarp/ProteinsProtAA.pdb --proteintraj $out_folder/output_schnarp/ProteinsProtAA.mdcrd --outputtop $out_folder/output_schnarp/Proteinspredicted.pdb --outputtraj $out_folder/output_schnarp/Proteinspredicted.mdcrd\n";
                `python $appFolder/Rebuild_nmerProteins.py --ref $out_folder/output_schnarp/reference.pdb --inputtop input/ProteinsChains.pdb --inputtraj $out_folder/output_schnarp/cgenarate_traj.mdcrd --proteintop $out_folder/output_schnarp/ProteinsProtAA.pdb --proteintraj $out_folder/output_schnarp/ProteinsProtAA.mdcrd --outputtop $out_folder/output_schnarp/Proteinspredicted.pdb --outputtraj $out_folder/output_schnarp/Proteinspredicted.mdcrd`;

                # Converting trajectory format for NGL
                #`mdconvert $out_folder/output_schnarp/cgenarate_traj.mdcrd -o $out_folder/output_schnarp/display/traj.dcd`;
                # Cpptraj: From mdcrd trajectory to DCD trajectory.
                #print "cpptraj -p $out_folder/output_schnarp/display/traj.pdb -y $out_folder/output_schnarp/cgenarate_traj.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log\n";
                #`cpptraj -p $out_folder/output_schnarp/display/traj.pdb -y $out_folder/output_schnarp/cgenarate_traj.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log 2>&1`;      
                print "cpptraj -p $out_folder/output_schnarp/Proteinspredicted.pdb  -y $out_folder/output_schnarp/Proteinspredicted.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log\n";
                `cpptraj -p $out_folder/output_schnarp/Proteinspredicted.pdb -y $out_folder/output_schnarp/Proteinspredicted.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log 2>&1`;   

                # Fix PDB file with consecutive residue numbering
                `perl $scriptsFolder/fixPDB.pl $out_folder/output_schnarp/Proteinspredicted.pdb > $out_folder/output_schnarp/Proteinspredicted_fixed.pdb`;

                # Generate dummy PARMTOP file from PDB
                `perl $scriptsFolder/prepTOP.pl $out_folder/output_schnarp/Proteinspredicted_fixed.pdb > $out_folder/output_schnarp/struc.prmtop`;

                # Copying PDB and PRMTOP files to display folder
                `cp $out_folder/output_schnarp/Proteinspredicted_fixed.pdb $out_folder/output_schnarp/display/traj.pdb`;   
                `cp $out_folder/output_schnarp/struc.prmtop $out_folder/output_schnarp/display/struc.prmtop`;   
        }
        else{
                # Adding proteins...
                print "python $appFolder/addProteinsToTraj.py --inputtraj $out_folder/output_schnarp/cgenarate_traj.mdcrd --outputtop $out_folder/output_schnarp/ProteinsProtAA.pdb --outputtraj $out_folder/output_schnarp/ProteinsProtAA.mdcrd \n";
                `python $appFolder/addProteinsToTraj.py --inputtraj $out_folder/output_schnarp/cgenarate_traj.mdcrd --outputtop $out_folder/output_schnarp/ProteinsProtAA.pdb --outputtraj $out_folder/output_schnarp/ProteinsProtAA.mdcrd`;

                # Converting trajectory format for NGL
                #`mdconvert $out_folder/output_schnarp/cgenarate_traj.mdcrd -o $out_folder/output_schnarp/display/traj.dcd`;
                # Cpptraj: From mdcrd trajectory to DCD trajectory.
                #print "cpptraj -p $out_folder/output_schnarp/display/traj.pdb -y $out_folder/output_schnarp/cgenarate_traj.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log\n";
                #`cpptraj -p $out_folder/output_schnarp/display/traj.pdb -y $out_folder/output_schnarp/cgenarate_traj.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log 2>&1`;      
                print "cpptraj -p $out_folder/output_schnarp/ProteinsProtAA.pdb  -y $out_folder/output_schnarp/ProteinsProtAA.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log\n";
                `cpptraj -p $out_folder/output_schnarp/ProteinsProtAA.pdb -y $out_folder/output_schnarp/ProteinsProtAA.mdcrd -x $out_folder/output_schnarp/display/traj.dcd > cpptraj_dcd.log 2>&1`;   

                # Fix PDB file with consecutive residue numbering
                `perl $scriptsFolder/fixPDB.pl $out_folder/output_schnarp/ProteinsProtAA.pdb > $out_folder/output_schnarp/ProteinsProtAA_fixed.pdb`;

                # Generate dummy PARMTOP file from PDB
                `perl $scriptsFolder/prepTOP.pl $out_folder/output_schnarp/ProteinsProtAA_fixed.pdb > $out_folder/output_schnarp/struc.prmtop`;

                # Copying PDB file to display folder
                `cp $out_folder/output_schnarp/ProteinsProtAA_fixed.pdb $out_folder/output_schnarp/display/traj.pdb`;   
                `cp $out_folder/output_schnarp/struc.prmtop $out_folder/output_schnarp/display/struc.prmtop`;   

        }


}
