#!/usr/bin/perl -w

use strict;

# Input Parameters.
my $long=@ARGV;
if ($long!=3){
    print "Usage: perl $0 <input.pdb> <input.top> <input.dcd>\n";
    exit(0);
}
my ($pdb,$top,$dcd)=@ARGV;

# NACCESS binary
#my $naccess = "/home/NACCESS";
my $naccess = "/app/Scripts/MCDNA/SASA/NACCESS";

# AVG Fiber SASA
my %avg_fiber;

# NACCESS ones
$avg_fiber{'DA'} = 171.86;
$avg_fiber{'DC'} = 176.41;
$avg_fiber{'DG'} = 175.72;
$avg_fiber{'DT'} = 185.78;

# Nucleic Bases
my %nucleic_codes = (
        "A" => 1,
        "T" => 1,
        "C" => 1,
        "G" => 1,
        "U" => 1,
        "DA" => 1,
        "DT" => 1,
        "DC" => 1,
        "DG" => 1,
        "DH" => 1,
        "DF" => 1,
        "RU" => 1,
        "RC" => 1,
        "RA" => 1,
        "RG" => 1,
        "A5" => 1,
        "A3" => 1,
        "C5" => 1,
        "C3" => 1,
        "T5" => 1,
        "T3" => 1,
        "G5" => 1,
        "G3" => 1,
        "U5" => 1,
        "U3" => 1,
        "DA5" => 1,
        "DA3" => 1,
        "DC5" => 1,
        "DC3" => 1,
        "DT5" => 1,
        "DT3" => 1,
        "DG5" => 1,
        "DG3" => 1,
        "RU5" => 1,
        "RU3" => 1,
        "RC5" => 1,
        "RC3" => 1,
        "RA5" => 1,
        "RA3" => 1,
        "RG5" => 1,
        "RG3" => 1
);

mkdir("PDB") if (! -s "PDB");

open PDBIN,">cpptraj.pdb.in";
print PDBIN "trajin $dcd\n";
print PDBIN "strip \@H*\n";
print PDBIN "trajout PDB/sasa.pdb multi\n";
print PDBIN "go\n";
close PDBIN;

#`/home/MCDNA/amber16/bin/cpptraj14 $top < cpptraj.pdb.in > cpptraj.pdb.out 2>&1`;
`cpptraj $top < cpptraj.pdb.in > cpptraj.pdb.out 2>&1`;

chdir("PDB");
foreach my $file (`ls --color=never *.pdb.*`){
	chomp($file);

	# sasa.pdb.27
	my @arr = split '\.',$file;
	my $newfile = "$arr[0]-$arr[2].$arr[1]";
	`cp $file $newfile`;

	# Running NACCESS process
	`$naccess/naccess $newfile`;

	my $rsafile="$arr[0]-$arr[2].rsa";

	#REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
	#REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
	#RES DG5 A   1   270.75 -99.9 270.75 -99.9   0.00 -99.9 143.03 -99.9 127.72 -99.9
	#RES DA  A   2   187.17 -99.9 170.77 -99.9  16.40 -99.9  84.93 -99.9 102.23 -99.9
	#RES DT  A   3   148.24 -99.9 132.71 -99.9  15.54 -99.9  86.41 -99.9  61.84 -99.9
	#RES DT  A   4    43.07 -99.9  38.94 -99.9   4.14 -99.9  35.62 -99.9   7.45 -99.9
	#RES DA  A   5   102.37 -99.9  86.96 -99.9  15.41 -99.9  64.47 -99.9  37.90 -99.9
	#RES DC  A   6   143.67 -99.9 125.36 -99.9  18.32 -99.9  64.49 -99.9  79.18 -99.9

	open OUT,">../sasa.$arr[2].dat";
	open FILE,"$rsafile";
	while (<FILE>){
		next if (!/^RES/);
		my $id = substr($_,4,3);
		$id=~s/ //g;
		my $num = substr($_,9,4);
		$num=~s/ //g;
		my $code = "$id-$num";
		my $abs = substr($_,14,8);
		#print "$code $abs\n" if ($nucleic_codes{$id});

		if ($nucleic_codes{$id}) {
			if (defined $avg_fiber{$id}){
				my $diff = $abs - $avg_fiber{$id};
				print OUT "$code $abs $avg_fiber{$id} $diff\n";
			}
			else{
				print OUT "$code $abs\n";
			}
		}
	}
	close FILE;
	close OUT;
}
chdir("..");
