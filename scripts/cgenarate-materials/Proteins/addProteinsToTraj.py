#!/usr/bin/env python
# coding: utf-8

# We need to know
# 
# 1. DNA sequence
# 1. TF position
# 1. TF length
# 1. TF pdb

import argparse
import numpy as np
from dnatraj import duplex as dx
import mdtraj as mdt

# Setup argument parser
parser = argparse.ArgumentParser(description="Process DNA and protein trajectories.")
parser.add_argument("--inputtraj", required=True, help="Input trajectory file path")
parser.add_argument("--outputtop", required=True, help="Output topology file path")
parser.add_argument("--outputtraj", required=True, help="Output trajectory file path")

args = parser.parse_args()

# Read the file content
with open('mcdna_prot.config', 'r') as file:
    lines = file.readlines()

# Process the lines and format the desired output
string = ""
for line in lines:
    parts = line.split()
    formatted_line = f"{parts[0]} pdb_inputs/{parts[1]}.pdb {parts[2]}"
    string += formatted_line + "\n"

# Remove the trailing newline character
string = string.strip()

print(string)

referencetop = 'input/Proteins.pdb'[:-4]+'Chains.pdb'

inputtraj = args.inputtraj
outputtop = args.outputtop
outputtraj = args.outputtraj

pdb_list = string.split()[1::3]

dyads = np.array(string.split()[::3], int)
dyads_len = np.array(string.split()[2::3], int)

mypdb = mdt.load(referencetop)
cg_traj = mdt.load(inputtraj, top=mypdb.top)

sel = mypdb.topology.select('name =~ "C1."')
c1pdb = mypdb.atom_slice(sel)
c1_traj = cg_traj.atom_slice(sel)

dx_c1_traj = dx.ComplementaryDuplex(c1_traj)

AApdb = [None]*len(dyads)

for i in range(0, len(dyads)):
    AApdb[i] = mdt.load(pdb_list[i], top=pdb_list[i])
    AApdb[i].xyz = np.tile(AApdb[i].xyz, (len(c1_traj), 1, 1))

    current_traj = dx_c1_traj[dyads[i]-1:dyads[i]+dyads_len[i]].traj

    selAA = AApdb[i].topology.select('name =~ "C1."')
    selCG = current_traj.topology.select('name =~ "C1."')

    for j in range(len(c1_traj)):
        current_frame = current_traj[j]
        current_frameAA = AApdb[i][0]
        current_frameAA.superpose(current_frame, atom_indices=selAA, ref_atom_indices=selCG)
        AApdb[i].xyz[j] = current_frameAA.xyz

    AApdb[i] = mdt.Trajectory(AApdb[i].xyz, AApdb[i].topology)
    sel = AApdb[i].topology.select('protein')
    AAProteinpdb = AApdb[i].atom_slice(sel)

    if i == 0:
        joined = c1_traj.stack(AAProteinpdb)
    else:
        joined = joined.stack(AAProteinpdb)

joined[0].save(outputtop)
joined.save(outputtraj)

