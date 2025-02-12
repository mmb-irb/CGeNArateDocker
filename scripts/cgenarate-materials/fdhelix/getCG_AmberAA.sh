#!/bin/bash

#first, compile: "bash compilefdhelix.sh"

#if sequence is short: "./fdhelix.exe abdna AAAAAAAAAATTTTTTTTTTATATATATATAAAAAAAAAATTTTTTTTTT | grep C1 > ATAT.pdb"
#if sequence is in file: "bash getfdhelix.sh seqfile.txt ATAT.pdb"

var=$(cat $1)

fdhelix.exe abdna $var > "structure_000000_AA_fdh.pdb"
grep "C1'" "structure_000000_AA_fdh.pdb" > "structure_000000_CG.pdb"

echo "source leaprc.DNA.bsc1
mol = loadpdb structure_000000_AA_fdh.pdb
savepdb mol structure_000000_AA_leap.pdb
quit" > leap_test.in

tleap -f leap_test.in > leap_test.out

echo "parm structure_000000_AA_leap.pdb
trajin structure_000000_AA_leap.pdb
strip @H=
trajout structure_000000_AA.pdb
run" > cpptraj_test.in

cpptraj -i cpptraj_test.in > cpptraj_test.out

#rm leap${2}.in leap${2}.out cpptraj${2}.in cpptraj${2}.out ${2}AA_fdh.pdb ${2}AA_leap.pdb

