#!/usr/bin/perl -w

# Fix PDB output files generated with ptraj, with a maximum of 99.999 residues and 999.999 atoms.
# Solve column problems. Pdb doesn't necessarily has to has atoms/residues consecutively numbered.

use strict;

# Input Parameters.
my $long=@ARGV;
if ($long!=1){
    print "Usage: perl $0 <pdb_to_fix.pdb>\n";
    exit(0);
}
my ($pdb)=@ARGV;

my $atn = 1;
my $resn = 1;
my $nresAnt = 1;

open PDB,"$pdb";
while(<PDB>){
	chomp;
	if(/^ATOM/ or /^HETATM/){

		# Parsing Res Number and Atom Number
		my $orig_resn = substr($_,22,6);
		#my $orig_atn = substr($_,5,7);

		$resn++ if ($orig_resn != $nresAnt);

		my $at_txt = sprintf("%5d",$atn);
		my $res_txt = sprintf("%4d",$resn);

		my $tag = substr($_,0,6);
		my $lineINI = substr($_,11,11);
		my $lineEND = substr($_,26);

		#print "$_\n";
		print "$tag$at_txt$lineINI$res_txt$lineEND\n";

		$atn++;

		$nresAnt = $orig_resn;
	}
	elsif(/^TER/){
		print "$_\n";
	}
}
close PDB;

