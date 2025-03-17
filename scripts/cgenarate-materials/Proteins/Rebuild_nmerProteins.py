#!/usr/bin/env python
# coding: utf-8

import argparse
from dnatraj import duplex as dx
import pickle
import mdtraj as mdt
import numpy as np

# Setup argument parser
parser = argparse.ArgumentParser(description="Process DNA and protein trajectories.")
parser.add_argument("--ref", required=True, help="Reference AA file path")
parser.add_argument("--inputtop", required=True, help="Input topology file path")
parser.add_argument("--inputtraj", required=True, help="Input trajectory file path")
parser.add_argument("--proteintop", required=True, help="Output topology file path")
parser.add_argument("--proteintraj", required=True, help="Output trajectory file path")
parser.add_argument("--outputtop", required=True, help="Output AA topology file path")
parser.add_argument("--outputtraj", required=True, help="Output AA trajectory file path")

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

#referencetop = 'input/ProteinsAA.pdb'

referencetop = args.ref
inputtop = args.inputtop
inputtraj = args.inputtraj
proteintop = args.proteintop
proteintraj = args.proteintraj
outputtop = args.outputtop
outputtraj = args.outputtraj

# bp=40
backbone_nmers=10

#bases = bp*2
#firsttestbp = bp-backbone_nmers-1
#trainnmers = firsttestbp-backbone_nmers

#trainbp = (firsttestbp-1)
#trainbases = 2*trainbp

#testbp = backbone_nmers
#testbases = 2*testbp


# # Rebuild full oligo
pickFolder = "/app/Scripts/cgenarate-materials/Proteins/"
with open(f"{pickFolder}/transformer_backbone.pickle", "rb") as input_file:
    transformer_backbone=pickle.load(input_file)
with open(f"{pickFolder}/transformerA.pickle", "rb") as input_file:
    transformerA=pickle.load(input_file)
with open(f"{pickFolder}/transformerC.pickle", "rb") as input_file:
    transformerC=pickle.load(input_file)

mypdb = mdt.load(referencetop)

sel = mypdb.topology.select('not type H')
mynoHpdb = mypdb.atom_slice(sel)

# mynoHpdb.save(outputtop)

sel = mypdb.topology.select('name =~ "C1."')
myc1pdb = mypdb.atom_slice(sel)

cg_test_backbone = mdt.load(inputtraj, top=inputtop)#[::10]
dx_test_backbone = dx.ComplementaryDuplex(cg_test_backbone)
lendx=len(dx_test_backbone)

dx_test_backbone.sequence#,dx_test_backbone.sequence==Sequences[simnum-1]

# cg_test_backbone.save("AtomAAAA_skip.mdcrd")

#Topology with extra phosphates for full reconstruction
auxtop=mynoHpdb.top.copy()

for name in ("OP2", "OP1", "P"):
    auxtop.insert_atom(name=name,
                       element=next(auxtop.atoms_by_name(name)).element,
                       residue=auxtop._residues[0],
                       index=0,
                       rindex=0)
    
for name in ("OP2", "OP1", "P"):
    auxtop.insert_atom(name=name,
                       element=next(auxtop.atoms_by_name(name)).element,
                       residue=auxtop._residues[lendx],
                       index=auxtop._residues[lendx]._atoms[0].index,
                       rindex=0) 

auxpdb=mdt.Trajectory(np.zeros((1,auxtop._numAtoms,3)),auxtop)

# ## Rebuild backbone

# ### Rebuild10mers

dx_nmer_pred_backbone = []

sel = auxpdb.topology.select('name =~ ".*\'.*|.*P.*"') 
dx_mypdb_backbone = dx.ComplementaryDuplex(auxpdb.atom_slice(sel))

for i in range(lendx // backbone_nmers):

    dx_backbone = dx_test_backbone[backbone_nmers * i:backbone_nmers * (i + 1)]
    
    tmp_mean = dx_backbone.traj.xyz.mean(axis=1).reshape(-1,1,3)
    xfg_predtest_backbone = transformer_backbone.transform(
        dx_backbone.traj.xyz -tmp_mean) +tmp_mean

    dx_mypdbs = dx_mypdb_backbone[backbone_nmers * i:backbone_nmers * (i + 1)]
    dx_pred = mdt.Trajectory(xfg_predtest_backbone, dx_mypdbs.traj.topology)
    dx_nmer_pred_backbone.append(dx.ComplementaryDuplex(dx_pred))
    print(i)

if (lendx % backbone_nmers != 0):

    dx_backbone = dx_test_backbone[lendx - backbone_nmers:lendx]

    tmp_mean = dx_backbone.traj.xyz.mean(axis=1).reshape(-1,1,3)
    xfg_predtest_backbone = transformer_backbone.transform(
        dx_backbone.traj.xyz -tmp_mean) +tmp_mean

    dx_mypdbs = dx_mypdb_backbone[lendx - backbone_nmers:lendx]
    dx_pred = mdt.Trajectory(xfg_predtest_backbone, dx_mypdbs.traj.topology)
    dx_nmer_pred_backbone.append(dx.ComplementaryDuplex(dx_pred))
    print("End")


# ### Join 10mers

cut = (backbone_nmers - lendx) % backbone_nmers
dx_last_pred_backbone = dx_nmer_pred_backbone[-1][cut:backbone_nmers]
dx_pred_backbone = dx.sstack(dx_nmer_pred_backbone[:-1] +
                             [dx_last_pred_backbone])
fg_predtest_backbone = dx_pred_backbone.traj

# fg_predtest_backbone.save_pdb('test_backbone.pdb')


# ## Rebuild Adenine from backbone 

# Rebuilt Adenine in test from "fg_predtest_backbone"

Aindexes_test_backbone = [
    res.index for res in cg_test_backbone.topology._residues
    if res.name[:2] == 'DA'
]

cg_test_Afrombackbone = [None] * len(Aindexes_test_backbone)

for i, index in enumerate(Aindexes_test_backbone):
    vmdsel = 'resid ' + str(index)
    sel = fg_predtest_backbone.topology.select(vmdsel)
    tempA = fg_predtest_backbone.atom_slice(sel)

    vmdsel = 'resid ' + str(lendx * 2 - 1 - index)
    sel = fg_predtest_backbone.topology.select(vmdsel)
    tempT = fg_predtest_backbone.atom_slice(sel)

    cg_test_Afrombackbone[i] = tempA.stack(tempT, True)

testtop=auxtop

fg_predtest_Afrombackbone = [None]*len(Aindexes_test_backbone)
print(len(Aindexes_test_backbone))

for i, index in enumerate(Aindexes_test_backbone):
    
    vmdsel = 'resid ' + str(index)
    sel = testtop.select(vmdsel)
    topA = testtop.subset(sel)

    vmdsel = 'resid ' + str(lendx*2-1-index)
    sel = testtop.select(vmdsel)
    topT = testtop.subset(sel)
    
    topAT = topA.join(topT)
    
    tmp_mean = cg_test_Afrombackbone[i].xyz.mean(axis=1).reshape(-1,1,3)
    xfg_predtest_Afrombackbone = transformerA.transform(cg_test_Afrombackbone[i].xyz -tmp_mean) +tmp_mean
    #merge with topology
    fg_predtest_Afrombackbone[i] = mdt.Trajectory(xfg_predtest_Afrombackbone, topAT) 
    
    
#     sel1 = fg_predtest_Afrombackbone[i].topology.select('name =~ ".*\'.*|.*P.*"')
#     sel2 = cg_test_Afrombackbone[i].topology.select('all')
#     for j in range(len(fg_predtest_Afrombackbone[i])):     
#         fg_predtest_Afrombackbone[i][j].superpose(cg_test_Afrombackbone[i][j],atom_indices=sel1,ref_atom_indices=sel2)
    # print(index)

# def mdtstack (trajectories, keep_resSeq=True):
#     #trajectories to stack (same time)
    
#     if len(trajectories) < 2:
#         print('why')
#         return
    
#     finaltraj=trajectories[0].stack(trajectories[1], keep_resSeq=True)
    
#     for i in trajectories[2:]:
#         finaltraj=finaltraj.stack(i, keep_resSeq=True)
        
#     return finaltraj

# mdtstack(fg_predtest_Afrombackbone, keep_resSeq=True).save_pdb('predtest_Afrombackbone.pdb')


# ## Rebuild Cytosine from backbone 

# Rebuilt Cytosine in test from "fg_predtest_backbone"

#select based on original traj
Cindexes_test_backbone = [res.index for res in cg_test_backbone.topology._residues if res.name[:2] == 'DC'] 

cg_test_Cfrombackbone = [None]*len(Cindexes_test_backbone)

for i, index in enumerate(Cindexes_test_backbone):
    vmdsel = 'resid ' + str(index)
    sel = fg_predtest_backbone.topology.select(vmdsel)
    tempC = fg_predtest_backbone.atom_slice(sel)
    
    vmdsel = 'resid ' + str(lendx*2-1-index)
    sel = fg_predtest_backbone.topology.select(vmdsel)
    tempG = fg_predtest_backbone.atom_slice(sel)
    
    cg_test_Cfrombackbone[i]=tempC.stack(tempG,True)

testtop=auxtop

fg_predtest_Cfrombackbone = [None]*len(Cindexes_test_backbone)
print(len(Cindexes_test_backbone))

for i, index in enumerate(Cindexes_test_backbone):
    
    vmdsel = 'resid ' + str(index)
    sel = testtop.select(vmdsel)
    topC = testtop.subset(sel)

    vmdsel = 'resid ' + str(lendx*2-1-index)
    sel = testtop.select(vmdsel)
    topG = testtop.subset(sel)
    
    topCG = topC.join(topG)
    
    tmp_mean = cg_test_Cfrombackbone[i].xyz.mean(axis=1).reshape(-1,1,3)
    xfg_predtest_Cfrombackbone = transformerC.transform(cg_test_Cfrombackbone[i].xyz -tmp_mean) +tmp_mean
    #merge with topology
    fg_predtest_Cfrombackbone[i] = mdt.Trajectory(xfg_predtest_Cfrombackbone, topCG)
    
#     sel1 = fg_predtest_Cfrombackbone[i].topology.select('name =~ ".*\'.*|.*P.*"')
#     sel2 = cg_test_Cfrombackbone[i].topology.select('all')
#     for j in range(len(fg_predtest_Cfrombackbone[i])):
#         fg_predtest_Cfrombackbone[i][j].superpose(cg_test_Cfrombackbone[i][j],atom_indices=sel1,ref_atom_indices=sel2)
    # print(index)

# def mdtstack (trajectories, keep_resSeq=True):
#     #trajectories to stack (same time)
    
#     if len(trajectories) < 2:
#         print('why')
#         return
    
#     finaltraj=trajectories[0].stack(trajectories[1], keep_resSeq=True)
    
#     for i in trajectories[2:]:
#         finaltraj=finaltraj.stack(i, keep_resSeq=True)
        
#     return finaltraj

# mdtstack(fg_predtest_Cfrombackbone, keep_resSeq=True).save_pdb('predtest_Cfrombackbone.pdb')

# ## Rejoin A and C predictions in backbone

tostack=[]
Astack=iter(fg_predtest_Afrombackbone)
Tstack=iter(fg_predtest_Afrombackbone[::-1])
Cstack=iter(fg_predtest_Cfrombackbone)
Gstack=iter(fg_predtest_Cfrombackbone[::-1])

for i in range(lendx):
    if i in Aindexes_test_backbone:
        tostack.append(dx.ComplementaryDuplex(next(Astack)))
    elif lendx*2-1-i in Aindexes_test_backbone:
        tostack.append(dx.ComplementaryDuplex(next(Tstack)).invert())
    elif i in Cindexes_test_backbone:
        tostack.append(dx.ComplementaryDuplex(next(Cstack)))
    elif lendx*2-1-i in Cindexes_test_backbone:
        tostack.append(dx.ComplementaryDuplex(next(Gstack)).invert())
    else:
        print("Error",i)

fg_predtest_0 = dx.sstack(tostack) 

# fg_predtest_0.traj.save_pdb('test_predicted0.pdb')

sel = fg_predtest_0.traj.topology.select(
    'not ((resid 0 %d) and (name =~ ".*P.*"))' % (lendx))
fg_predtest = fg_predtest_0.traj.atom_slice(sel)

print("Backmapping complete")

pdb_list = string.split()[1::3]

dyads = np.array(string.split()[::3],int)
dyads_len = np.array(string.split()[2::3],int)

mypdb = mdt.load(inputtop)
cg_traj = mdt.load(inputtraj, top=mypdb.top)

sel = mypdb.topology.select('name =~ "C1."')
c1pdb = mypdb.atom_slice(sel)
c1_traj = cg_traj.atom_slice(sel)

dx_c1_traj = dx.ComplementaryDuplex(c1_traj)
AApdb = [None]*len(dyads)

for i in range(0,len(dyads)):
    AApdb[i] = mdt.load(pdb_list[i], top=pdb_list[i])
    AApdb[i].xyz = np.tile(AApdb[i].xyz, (len(c1_traj), 1, 1))

    current_traj = dx_c1_traj[dyads[i]-1:dyads[i]+dyads_len[i]].traj

    selAA = AApdb[i].topology.select('name =~ "C1."')
    selCG = current_traj.topology.select('name =~ "C1."')

    for j in range(len(c1_traj)):
        current_frame = current_traj[j]
        current_frameAA = AApdb[i][0]
        current_frameAA.superpose(current_frame,atom_indices=selAA,ref_atom_indices=selCG)
        # AApdb[i] = AApdb[i].join(current_frameAA, check_topology=False)

        AApdb[i].xyz[j] = current_frameAA.xyz
        # print(j)

    AApdb[i] = mdt.Trajectory(AApdb[i].xyz,AApdb[i].topology)
    # AApdb[i].time = np.arange(len(AApdb[i].xyz))
    # print(AApdb[i].time,len(AApdb[i].xyz))
    sel = AApdb[i].topology.select('protein')
    AAProteinpdb = AApdb[i].atom_slice(sel)

    if i == 0:
        # joined = c1_traj.stack(AAProteinpdb)
        joined = AAProteinpdb
    else:
        joined = joined.stack(AAProteinpdb)
    print(i,pdb_list[i])

print('Protein complete')
# joined[0].save(outputtop)
# joined.save(outputtraj)

# proteinpdb = mdt.load(proteintop)
# protein_traj = mdt.load(proteintraj, top=proteinpdb.top)

Full = fg_predtest.stack(joined)
# Full = fg_predtest

Full.save(outputtraj)

Full[0].save(outputtop)

print('Printing complete')
