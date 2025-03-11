#!/usr/bin/perl -w

use strict;

# Script Input Parameters
my $long=@ARGV;
if ($long < 6 or $long > 7){
    print "Usage: perl $0 <seq_file> <num_structures> <segment_length> <linking_number> <number_iterations> <all-atom> [0 (or
 None): <eqStruct>, 1: <create_traj>, 2: <create_both_eqStruct_and_traj>]\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 5 -1 5000 1 0  # To generate a single structure (AA) (equilibrated)\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 5 -1 5000 1 1  # To generate an all-atom trajectory\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 5 -1 5000 0 1  # To generate a Coarse-Grained trajectory\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 5 -1 5000 1 2  # To generate both: AA structure + trajectory\n";
    print "Example: perl $0 seq_56merSPCE.dat 10 5 -1 5000 0 2  # To generate both: CG structure + trajectory\n";
    exit(0);
}

# Input Parameters.
my ($seq,$nstructs,$seglen,$deltalk,$niterations,$rebuilt,$traj)=@ARGV;

$traj = 0 if !(defined $traj);

my $system = "STRUCT";
$system = "TRAJ" if ($traj == 1);
$system = "STRUCT_AND_TRAJ" if ($traj == 2);

my $resolution = "AA";
$resolution = "CG" if (!$rebuilt);

# Auxiliar scripts/binaries

# MCDNA binaries folder. Execution must be done inside the folder, 
# as it needs auxiliar files hardcoded in the C code from this folder.
my $MCDNAbinStruct = "getfdhelix.sh";
my $MCDNAbinStruct_AA = "getCG_AmberAA.sh";

# CGenerate Paths.
my $MCDNAbinTraj = "/app/Scripts/cgenarate-Mitochondria";
my $MCDNAfdhelix = "/app/Scripts/cgenarate-materials/fdhelix";

# GLIMPS Path.
my $glimps = "/app/Scripts/cgenarate-materials/GLIMPS";

# MCDNA 2025 extended PATH.
$ENV{'PATH'}="$ENV{'PATH'}:/app/Scripts/cgenarate-materials/fdhelix/";

# CONDA extended PATH.
$ENV{'PATH'}="$ENV{'PATH'}:/opt/conda/bin";

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

# Sequence length *2 (2 strands)
my $num = length($arr);
$num *= 2;

# Initializing output folders.
my $curr_folder = `pwd`;
chomp($curr_folder);

if($traj == 1 or $traj == 2){

	# GENERATING TRAJECTORY (ENSEMBLE with $nstructs STRUCTURES)
        print "## GENERATING TRAJECTORY (ENSEMBLE with $nstructs STRUCTURES)\n";

	my $out_folder = "$curr_folder/TRAJ_$resolution";
	#mkdir("$out_folder") if (! -s "$out_folder");
	chdir("$out_folder");

	# STEP 1: Executing fdhelix (generating structure from sequence)
	print "# STEP 1: Executing fdhelix...\n";
	&runMCDNA_new($MCDNAbinStruct,$out_folder,$seq,$nstructs,$rebuilt);

	# STEP 2: Executing CGeNArate
	# It needs to be executed in the same folder where the binary is placed (auxiliar files hardcoded)
	print "# STEP 2: Executing CGenerate...\n";
	&runCgenerate($MCDNAbinTraj,$out_folder,$seq,$nstructs,$rebuilt);

	print "# STEP 3: Building trajectory...\n";
	
	mkdir("$out_folder/output_schnarp") if (! -s "$out_folder/output_schnarp");
	mkdir("$out_folder/output_schnarp/display") if (! -s "$out_folder/output_schnarp/display");
	chdir("$out_folder/output_schnarp");

	open PTRAJIN,">reference_pdb_cpptraj.in";
	#print PTRAJIN "trajin $out_folder/output_pdb/out.mdcrd 1 1 1\n";
	print PTRAJIN "trajin $out_folder/output_schnarp/cgenarate_traj.mdcrd 1 1 1\n";
	print PTRAJIN "trajout mc_dna_str.pdb pdb \n";
	print PTRAJIN "go\n";
	print PTRAJIN "quit\n";
	close PTRAJIN;

	open PTRAJIN,">make_traj_cpptraj.in";
	print PTRAJIN "trajin $out_folder/output_schnarp/cgenarate_traj.mdcrd\n";
	print PTRAJIN "reference $out_folder/output_schnarp/str.pdb\n";
	print PTRAJIN "rms reference :1-$num\n";
	print PTRAJIN "trajout cgenarate_traj.fitted.mdcrd \n";
	print PTRAJIN "go\n";
	print PTRAJIN "quit\n";
	close PTRAJIN;

	if (! -s "$out_folder/output_pdb/structure_000000_CG.pdb") {
		`cp $out_folder/output_pdb/structure_000000_circ.pdb $out_folder/output_pdb/structure_000000_CG.pdb`;
	}

		# Cpptraj: Generate reference PDB.
	print "cpptraj -p $out_folder/output_pdb/structure_000000_CG.pdb < reference_pdb_cpptraj.in > reference_pdb_cpptraj.log\n";
	`cpptraj -p $out_folder/output_pdb/structure_000000_CG.pdb < reference_pdb_cpptraj.in > reference_pdb_cpptraj.log 2>&1`;

		# Cpptraj: Generate RMS fitted mdcrd trajectory.
	print "cpptraj -p $out_folder/output_pdb/structure_000000_CG.pdb < make_traj_cpptraj.in > cpptraj.log\n";
	`cpptraj -p $out_folder/output_pdb/structure_000000_CG.pdb < make_traj_cpptraj.in > cpptraj.log 2>&1`;

	# Cpptraj: From mdcrd trajectory to DCD trajectory.
	print "cpptraj -p $out_folder/output_pdb/structure_000000_CG.pdb -y $out_folder/output_schnarp/cgenarate_traj.fitted.mdcrd -x $out_folder/output_schnarp/display/mc_dna_str.dcd > cpptraj_dcd.log\n";
	`cpptraj -p $out_folder/output_pdb/structure_000000_CG.pdb -y $out_folder/output_schnarp/cgenarate_traj.fitted.mdcrd -x $out_folder/output_schnarp/display/mc_dna_str.dcd > cpptraj_dcd.log 2>&1`;

	if ($resolution eq "AA"){
	
		print "# STEP 4: Atomistic Reconstruction (GLIMPS)...\n";

                # GLIMPS atomistic reconstruction
       	        print "python $glimps/Rebuild_nmer_Circ_web.py $out_folder/output_schnarp/str.pdb $out_folder/output_schnarp/display/mc_dna_str.dcd $out_folder/output_schnarp/display/traj_glimps.pdb $out_folder/output_schnarp/display/traj_glimps.dcd\n";
               	`python $glimps/Rebuild_nmer_Circ_web.py $out_folder/output_schnarp/str.pdb $out_folder/output_schnarp/display/mc_dna_str.dcd $out_folder/output_schnarp/display/traj_glimps.pdb $out_folder/output_schnarp/display/traj_glimps.dcd`;
			
		`cp $out_folder/output_schnarp/display/traj_glimps.pdb $out_folder/output_schnarp/display/traj.pdb`;
		`cp $out_folder/output_schnarp/display/traj_glimps.dcd $out_folder/output_schnarp/display/traj.dcd`;

		# Converting trajectory format for Curves+
		`mdconvert $out_folder/output_schnarp/display/traj.dcd -o $out_folder/output_schnarp/display/traj.nc`;
	}
	else{
		`cp $out_folder/output_schnarp/str.pdb $out_folder/output_schnarp/display/traj.pdb`;
		`cp $out_folder/output_schnarp/display/mc_dna_str.dcd $out_folder/output_schnarp/display/traj.dcd`;
	}
}

if($traj == 2 or $traj == 0){

	# GENERATING STRUCTURE 
        print "## GENERATING STRUCTURE \n";

	my $out_folder = "$curr_folder/EQ_$resolution";
	chdir("$out_folder");

	# STEP 1: Executing CGeNArate
	# It needs to be executed in the same folder where the binary is placed (auxiliar files hardcoded)
	print "# STEP 1: Executing CGeNArate...\n";
	&runMCDNA_new($MCDNAbinStruct,$out_folder,$seq,$nstructs,$rebuilt);

	if ($resolution eq "AA"){

		# AA STRUCTURE

		print "# STEP 3: Building trajectory...\n";
	
		mkdir("$out_folder/output_schnarp") if (! -s "$out_folder/output_schnarp");
		chdir("$out_folder/output_schnarp");

		#open PTRAJIN,">reference_pdb_cpptraj.in";
		#print PTRAJIN "trajin $out_folder/output_pdb/out.mdcrd 1 1 1\n";
		#print PTRAJIN "trajout mc_dna_str.pdb pdb \n";
		#print PTRAJIN "go\n";
		#print PTRAJIN "quit\n";
		#close PTRAJIN;

		# Cpptraj: Generate reference PDB.
		#print ". /opt/conda/etc/profile.d/conda.sh; conda activate glimps_env; cpptraj -p $out_folder/output_pdb/out.top < reference_pdb_cpptraj.in > reference_pdb_cpptraj.log\n";
		#`. /opt/conda/etc/profile.d/conda.sh; conda activate glimps_env; cpptraj -p $out_folder/output_pdb/out.top < reference_pdb_cpptraj.in > reference_pdb_cpptraj.log 2>&1`;

		# Standard output file name (for the web server)
		#`cp mc_dna_str.pdb str.pdb`;
		#`cp $out_folder/output_pdb/structure_000000_AA.pdb str.pdb`;
		`cp $out_folder/output_pdb/structure_000000_AA_fdh_circ.pdb str.pdb`;
	}
}

# Building downloadable tar package
print "Building downloadable tar package...\n";
chdir("$curr_folder");
#`tar -zcvf mcdna.tgz EQ*` if ($traj == 0);
#`tar -zcvf mcdna.tgz TRAJ*` if ($traj == 1);
#`tar -zcvf mcdna.tgz EQ* TRAJ*` if ($traj == 2);

my $method = "mcdna";
&download_readme($arr,$method,$resolution,$system);
if($traj == 1 or $traj == 2){
	&download_traj($resolution);
}
elsif($traj == 0 or $traj == 2){
	&download_eq($resolution);
}

chdir("$curr_folder/download");
my $tf = "$method-$resolution-$system".".tgz";
`tar -zcvf $tf *`;
chdir("..");
`mv download/$tf .`;

# Cleaning downloadable (and duplicated) info
#`rm -r download`;

print "DONE!!\n";

##################################################################################################################

sub runCgenerate {

	my ($dir_binary,$out_folder,$seq,$nstructs,$atom) = @_;

	mkdir("$out_folder") if (!-s "$out_folder");
	mkdir("$out_folder/output_pdb") if (!-s "$out_folder/output_pdb");
	mkdir("$out_folder/output_schnarp") if (!-s "$out_folder/output_schnarp");

	# CGeNArate only works if executed where the binary is (uses hard-coded auxiliary files)
	chdir("$dir_binary");

	my $config = "$out_folder/output_schnarp/cgenerate_config.toml";
        open CONFIG,">$config";
        print CONFIG "# CGeNArate config file, created by CGeNArate web server #\n";
	print CONFIG "[simulation]\n";
	print CONFIG "name = 'CGeNArate web server - Circular DNA'\n";
	print CONFIG "T0 = 298.0 # Temperature (K)\n";
	print CONFIG "dt0 = 0.1 # timestep in picoseconds\n";
	print CONFIG "N = $num # Number of beads\n";
	print CONFIG "frames = $nstructs # Number of frames\n";
	print CONFIG "steps = 200 # Number of steps between frames\n";
	print CONFIG "linear = false\n";
	print CONFIG "seed = 1234 # random seed\n\n";
	print CONFIG "[input]\n";
	print CONFIG "file = '$out_folder/output_schnarp/str.pdb'\n\n";
	print CONFIG "[output]\n";
	print CONFIG "coordinates = '$out_folder/output_schnarp/cgenarate_traj.mdcrd'\n";
	print CONFIG "energy = '$out_folder/output_schnarp/cgenarate_energy.csv'\n";
	print CONFIG "plots = '$out_folder/output_schnarp/cgenarate_plots.csv'\n\n";
	print CONFIG "[restore]\n";
	print CONFIG "mode = 'print'\n";
	print CONFIG "current = '$out_folder/output_schnarp/current.pdb'\n";
	print CONFIG "restart = '$out_folder/output_schnarp/restart.pdb'\n";
	close CONFIG;

	print "$dir_binary/CGeNArate.seq.exe $config\n";
	`$dir_binary/CGeNArate.seq.exe $config`;
	
	if($atom){
		# Run GLIMPS
	}
	else{
	}
}

sub runMCDNA_new {

	my ($dir_binary,$out_folder,$seq,$nstructs,$atom) = @_;

	mkdir("$out_folder") if (!-s "$out_folder");
	mkdir("$out_folder/output_pdb") if (!-s "$out_folder/output_pdb");
	mkdir("$out_folder/output_schnarp") if (!-s "$out_folder/output_schnarp");

	chdir("$out_folder/output_pdb");
	if($atom){
		print "$MCDNAbinStruct_AA $seq $out_folder/output_pdb/structure_000000_ \n";
		`$MCDNAbinStruct_AA $seq ${out_folder}/output_pdb/structure_000000_`;

		print "python $MCDNAfdhelix/buildCircular.py $out_folder/output_pdb/structure_000000_AA_fdh.pdb $out_folder/output_pdb/structure_000000_AA_fdh_circ.pdb $deltalk\n"; 
		`python $MCDNAfdhelix/buildCircular.py $out_folder/output_pdb/structure_000000_AA_fdh.pdb $out_folder/output_pdb/structure_000000_AA_fdh_circ.pdb $deltalk`;

                # Standard output file name (for the web server)
		`cp $out_folder/output_pdb/structure_000000_AA_fdh_circ.pdb $out_folder/output_pdb/str.pdb`;
		`cp $out_folder/output_pdb/structure_000000_AA_fdh_circ.pdb $out_folder/output_schnarp/str.pdb`;
	}
	else{
		print "$MCDNAbinStruct $seq $out_folder/output_pdb/structure_000000.pdb \n";
		`$MCDNAbinStruct $seq $out_folder/output_pdb/structure_000000.pdb`;

		print "python $MCDNAfdhelix/buildCircular.py $out_folder/output_pdb/structure_000000.pdb $out_folder/output_pdb/structure_000000_circ.pdb $deltalk\n"; 
		`python $MCDNAfdhelix/buildCircular.py $out_folder/output_pdb/structure_000000.pdb $out_folder/output_pdb/structure_000000_circ.pdb $deltalk`;

                # Standard output file name (for the web server)
                `cp $out_folder/output_pdb/structure_000000_circ.pdb $out_folder/output_schnarp/str.pdb`;
	}
}

sub download_readme {
        my ($seq,$method,$resolution,$system) = @_;

        mkdir("download") if (! -s "download");

        #my $nstructs = 0;
        #my $log = "TRAJ_$resolution/output_schnarp/display/cpptraj_dcd.log";
        #if (-s "$log"){
        #        my $g = `grep "frames and processed" $log`;
        #        $g =~ /Read (\d+) frames and processed/;
        #        $nstructs = $1;
        #}

        # Auxiliary Files
        open SUM,">download/summary.txt";
        print SUM "# CGeNArate process summary #\n";
        print SUM "Sequence: $seq\n";
        print SUM "Method: $method\n";
        print SUM "Resolution: $resolution\n";
        print SUM "System: $system\n";
        #print SUM "Number of structures: $nstructs\n" if ($system =~ /TRAJ/ and $nstructs);
        print SUM "Number of structures: $nstructs\n"; 
        close SUM;
}

sub download_eq {
        my ($resolution) = @_;

        mkdir("download") if (! -s "download");

        # MCDNA Structure
        mkdir("download/STRUCT") if (! -s "download/STRUCT");
        mkdir("download/STRUCT/CGeNArate") if (! -s "download/STRUCT/CGeNArate");

        `cp EQ_$resolution/output_schnarp/str.pdb download/STRUCT/CGeNArate`;
        #`cp -r EQ_$resolution/output_helpar download/STRUCT/CGeNArate`;
	#`mv download/STRUCT/CGeNArate/output_helpar download/STRUCT/CGeNArate/bps_parms`;
}

sub download_traj {
        my ($resolution) = @_;

        mkdir("download") if (! -s "download");

        # CGeNArate Structure
        mkdir("download/TRAJ") if (! -s "download/TRAJ");
        mkdir("download/TRAJ/CGeNArate") if (! -s "download/TRAJ/CGeNArate");

        `cp TRAJ_$resolution/output_schnarp/display/traj.pdb download/TRAJ/CGeNArate`;
        `cp TRAJ_$resolution/output_schnarp/display/traj.dcd download/TRAJ/CGeNArate`;
        #`cp -r TRAJ_$resolution/output_helpar download/TRAJ/CGeNArate`;
	#`mv download/TRAJ/CGeNArate/output_helpar download/TRAJ/CGeNArate/bps_parms`;
}

