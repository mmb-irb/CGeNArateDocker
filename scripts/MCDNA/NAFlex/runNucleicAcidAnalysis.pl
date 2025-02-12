#!/usr/bin/perl 

use lib "/app/Scripts/MCDNA/NAFlex";
use NA_MDWeb_Common;
#use Switch;
use Switch 'fallthrough';
use strict;

require "NA_MDWeb_curves.pl";
require "NA_MDWeb_stiffness.pl";
require "NA_MDWeb_pcazip.pl";
require "NA_MDWeb_DistContactMapsProts.pl";

my $long=@ARGV;
if ($long < 6){
    #print "Usage: perl $0 <pdbFile> <topFile> <crdFile> <workdir> <output> <operation> [force] [netcdf]\n";
    print "Usage: perl $0 <pdbFile> <topFile> <crdFile> <workdir> <output> <operation/s>\n";
    print "Operations:\n";
    print "\tALL\n\tCurves\n\tStiffness\n\tPcazip\n\tDistanceContactMaps\n";
    exit(0);
}
my ($mypdb, $mytop, $mycrd, $workdir, $prefix, @operation) =@ARGV;

# AMBER folder and Environmental Vars.
#my $achome = "/home/MCDNA/amber16";
#$ENV{'AMBERHOME'} = "/home/MCDNA/amber16";

# setenv LD_LIBRARY_PATH 
#$ENV{'LD_LIBRARY_PATH'} = "$achome/lib:/opt/software/lib:/usr/lib/atlas:$ENV{'LD_LIBRARY_PATH'}";
#print "$ENV{'LD_LIBRARY_PATH'}\n";

# setenv PATH 
#$ENV{'PATH'} = "$achome/bin:/opt/software/bin:$ENV{'PATH'}";
#print "$ENV{'PATH'}\n";

my $scriptsdir = "/app/Scripts/MCDNA/NAFlex";
print "Workdir: $workdir\n";

my $prefix_orig = $prefix;

$prefix = 'MCDNA_prova' if (!$prefix);

my $pdb = `readlink -f $mypdb`;
chomp($pdb);
my $top = `readlink -f $mytop`;
chomp($top);
my $crd = `readlink -f $mycrd`;
chomp($crd);
my $crd_orig = $crd;
my $top_orig = $top;
my $pdb_orig = $pdb;

my $pdbNoH = $pdb;
my $crdNoH = $crd;

mkdir("$workdir/$prefix") if(! -e "$workdir/$prefix");
chdir("$workdir/$prefix");

my $NAFleWorkDir = "$workdir/$prefix";

# Generating completely stripped trajectory (no water nor ions)
mkdir("$NAFleWorkDir/INFO") if(! -e "$NAFleWorkDir/INFO");

print "Generating stripped traj...\n";
open IN,">$NAFleWorkDir/INFO/cpptraj.stripTrj.in";
print IN "trajin $crd\n";
print IN "strip :WAT\n";
print IN "strip :Cl-\n";
print IN "strip :Cl\n";
print IN "strip :Na+\n";
print IN "strip :K+\n";
print IN "strip :Mn++\n";
print IN "strip :MNG\n";
print IN "strip :FED\n";
print IN "strip :ALA,ARG,ASP,GLU,PHE,TYR,TRP,LYS,GLY,HIS,HIE,HID,HIP,CYS,CYX,SER,PRO,THR,MET,LEU,ILE,ASN,GLN,VAL\n";
print IN "autoimage\n";
print IN "rmsd !\@H* first out $NAFleWorkDir/INFO/prova.rmsd\n";
print IN "trajout $NAFleWorkDir/INFO/structure.stripped.trj\n";
close IN;
`cpptraj $top < $NAFleWorkDir/INFO/cpptraj.stripTrj.in > $NAFleWorkDir/INFO/cpptraj.stripTrj.out 2>&1`; 

print "Generating stripped traj without hydrogens (for pcazip)...\n";
open IN,">$NAFleWorkDir/INFO/cpptraj.stripTrj.noH.in";
print IN "trajin $crd\n";
print IN "strip :WAT\n";
print IN "strip :Cl-\n";
print IN "strip :Cl\n";
print IN "strip :Na+\n";
print IN "strip :K+\n";
print IN "strip :Mn++\n";
print IN "strip :MNG\n";
print IN "strip :FED\n";
print IN "strip :ALA,ARG,ASP,GLU,PHE,TYR,TRP,LYS,GLY,HIS,HIE,HID,HIP,CYS,CYX,SER,PRO,THR,MET,LEU,ILE,ASN,GLN,VAL\n";
print IN "strip \@H*\n";
print IN "autoimage\n";
print IN "trajout $NAFleWorkDir/INFO/structure.stripped.noH.trj\n";
close IN;
`cpptraj $top < $NAFleWorkDir/INFO/cpptraj.stripTrj.noH.in > $NAFleWorkDir/INFO/cpptraj.stripTrj.noH.out 2>&1`; 

print "Generating stripped top...\n";
open IN,">$NAFleWorkDir/INFO/cpptraj.stripTop.in";
print IN "parmstrip :WAT\n";
print IN "parmstrip :Cl-\n";
print IN "parmstrip :Cl\n";
print IN "parmstrip :Na+\n";
print IN "parmstrip :K+\n";
print IN "parmstrip :Mn++\n";
print IN "parmstrip :MNG\n";
print IN "parmstrip :FED\n";
print IN "parmstrip :ALA,ARG,ASP,GLU,PHE,TYR,TRP,LYS,GLY,HIS,HIE,HID,HIP,CYS,CYX,SER,PRO,THR,MET,LEU,ILE,ASN,GLN,VAL\n";
#print IN "parmstrip :1,20\n";
print IN "parmwrite out $NAFleWorkDir/INFO/structure.stripped.top amber\n";
print IN "go\n";
close IN;
#`cpptraj $top < $NAFleWorkDir/INFO/cpptraj.stripTop.in > $NAFleWorkDir/INFO/cpptraj.stripTop.out 2>&1`; 
`cp $top $NAFleWorkDir/INFO/structure.stripped.top`;

print "Generating stripped top without hydrogens (for pcazip) ...\n";
open IN,">$NAFleWorkDir/INFO/cpptraj.stripTop.noH.in";
print IN "parmstrip :WAT\n";
print IN "parmstrip :Cl-\n";
print IN "parmstrip :Cl\n";
print IN "parmstrip :Na+\n";
print IN "parmstrip :K+\n";
print IN "parmstrip :Mn++\n";
print IN "parmstrip :MNG\n";
print IN "parmstrip :FED\n";
print IN "parmstrip :ALA,ARG,ASP,GLU,PHE,TYR,TRP,LYS,GLY,HIS,HIE,HID,HIP,CYS,CYX,SER,PRO,THR,MET,LEU,ILE,ASN,GLN,VAL\n";
print IN "parmstrip \@H*\n";
print IN "parmwrite out $NAFleWorkDir/INFO/structure.stripped.noH.top amber\n";
print IN "go\n";
close IN;
#`cpptraj $top < $NAFleWorkDir/INFO/cpptraj.stripTop.noH.in > $NAFleWorkDir/INFO/cpptraj.stripTop.noH.out 2>&1`; 
`cp $top $NAFleWorkDir/INFO/structure.stripped.noH.top`;

print "Generating stripped pdb...\n";
open IN,">$NAFleWorkDir/INFO/cpptraj.stripPdb.in";
print IN "trajin $NAFleWorkDir/INFO/structure.stripped.trj 1 1 1\n";
print IN "trajout $NAFleWorkDir/INFO/structure.stripped.pdb\n";
close IN;
`cpptraj $NAFleWorkDir/INFO/structure.stripped.top < $NAFleWorkDir/INFO/cpptraj.stripPdb.in > $NAFleWorkDir/INFO/cpptraj.stripPdb.out 2>&1`; 

print "Generating stripped pdb without hydrogens (for pcazip)...\n";
open IN,">$NAFleWorkDir/INFO/cpptraj.stripPdb.noH.in";
print IN "trajin $NAFleWorkDir/INFO/structure.stripped.noH.trj 1 1 1\n";
print IN "trajout $NAFleWorkDir/INFO/structure.stripped.noH.pdb\n";
close IN;
`cpptraj $NAFleWorkDir/INFO/structure.stripped.noH.top < $NAFleWorkDir/INFO/cpptraj.stripPdb.noH.in > $NAFleWorkDir/INFO/cpptraj.stripPdb.noH.out 2>&1`; 

# Fixing possible problems with generated PDB:
my $gnum = `grep TER $NAFleWorkDir/INFO/structure.stripped.pdb | wc -l`;
$gnum +=0;
if ($gnum > 100){
	print "Problems with generated pdb ($gnum TERs), stripping TER tags...\n";
	`cp $NAFleWorkDir/INFO/structure.stripped.pdb $NAFleWorkDir/INFO/structure.stripped.pdb.tmp`;
	`grep -v TER $NAFleWorkDir/INFO/structure.stripped.pdb.tmp > $NAFleWorkDir/INFO/structure.stripped.pdb`;
}

# Fixing possible problems with generated PDB without hydrogens:
$gnum = `grep TER $NAFleWorkDir/INFO/structure.stripped.noH.pdb | wc -l`;
$gnum +=0;
if ($gnum > 100){
	print "Problems with generated pdb ($gnum TERs), stripping TER tags...\n";
	`cp $NAFleWorkDir/INFO/structure.stripped.noH.pdb $NAFleWorkDir/INFO/structure.stripped.noH.pdb.tmp`;
	`grep -v TER $NAFleWorkDir/INFO/structure.stripped.pdb.tmp > $NAFleWorkDir/INFO/structure.stripped.noH.pdb`;
}

$crd = "$NAFleWorkDir/INFO/structure.stripped.trj" if (-s "$NAFleWorkDir/INFO/structure.stripped.trj");
$top = "$NAFleWorkDir/INFO/structure.stripped.top" if (-s "$NAFleWorkDir/INFO/structure.stripped.top");
$pdb = "$NAFleWorkDir/INFO/structure.stripped.pdb" if (-s "$NAFleWorkDir/INFO/structure.stripped.pdb");
$pdbNoH = "$NAFleWorkDir/INFO/structure.stripped.noH.pdb"; 
$crdNoH = "$NAFleWorkDir/INFO/structure.stripped.noH.trj";

if (! -s "$NAFleWorkDir/INFO/structure.stripped.trj"){
	print "Backup plan... converting DCD to TRJ...\n";
	open IN,">$NAFleWorkDir/INFO/cpptraj.dcdtocrd.in";
	print IN "trajin $crd\n";
	print IN "trajout $NAFleWorkDir/INFO/structure.stripped.trj mdcrd\n";
	close IN;
	`cpptraj $top < $NAFleWorkDir/INFO/cpptraj.dcdtocrd.in > $NAFleWorkDir/INFO/cpptraj.dcdtocrd.out 2>&1`; 

	$crd = "$NAFleWorkDir/INFO/structure.stripped.trj" if (-s "$NAFleWorkDir/INFO/structure.stripped.trj");
}

foreach my $operation (@operation) {
    switch ($operation) {

        case /ALL|Curves/ {

		# case "Curves"
		print "Executing Curves Analysis (using Curves+ program).\n\n";
		print "run_curves($NAFleWorkDir,$prefix,$pdb_orig,$top,$crd)\n";
		run_curves($NAFleWorkDir,$prefix,$pdb_orig,$top,$crd);

		last if $operation eq 'Curves';
		
	}

        case /ALL|Stiffness/ {

		# case "Stiffness"
		print "Executing Stiffness Constants Analysis (using Curves+ program).\n\n";
		#print "run_stiffness($NAFleWorkDir,$prefix,$pdb,$top,$crd)\n";
		run_stiffness($NAFleWorkDir,$prefix,$pdb,$top,$crd);

		last if ($operation eq 'Stiffness');
	}

        case /ALL|Pcazip/ {

		# case "Pcazip"
		print "Executing Principal Component Analysis with pcazip program.\n\n";
		#print "run_pcazip($NAFleWorkDir,$prefix,$pdb,$crd)\n";
		run_pcazip($NAFleWorkDir,$prefix,$pdbNoH,$crdNoH);

		last if ($operation eq 'Pcazip');
	}

        case /ALL|DistanceContactMaps/ {

		# case "DistanceContactMaps"
		print "Executing Distance Contact-Maps Analysis.\n\n";
		#print "run_distContactMaps($NAFleWorkDir,$prefix,$pdb,$crd)\n";
		run_distContactMaps($NAFleWorkDir,$prefix,$pdb_orig,$top_orig,$crd_orig);

		last if ($operation eq 'DistanceContactMaps');
	}
    }
}
