#!/usr/bin/env python
# coding: utf-8

import sys
import dnatraj as dnat
import pickle
import mdtraj as mdt
import numpy as np
from mdplus.multiscale import Glimps
from dnatraj import duplex as dx

if len(sys.argv) != 5: 
	print("Usage: sys.argv[0] <reference_top> <input_traj> <output_top> <output_traj>") 
	sys.exit(1) 

referencetop = sys.argv[1] 
inputtraj = sys.argv[2] 
outputtop = sys.argv[3] 
outputtraj = sys.argv[4] 

#referencetop = 'input/' + names[i] + 'AA.pdb'
#inputtraj = 'input/' + names[i] + 'aligned.mdcrd'
#outputtop = 'output/' + names[i] + 'predicted.pdb'
#outputtraj = 'output/' + names[i] + 'predicted.mdcrd'

names = ["2lef", "AGCG", "AGCT", "CGTG", "1zgw", "1j5n", "56merL", "CTAG_flex"]
i = 6
i, names[i]

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

#dirGlimps = "/orozco/services/MCDNA_dev/Scripts/2025/cgenarate-materials/GLIMPS/"
dirGlimps = "/app/Scripts/cgenarate-materials/GLIMPS/"

with open(dirGlimps + "transformer_backbone.pickle", "rb") as input_file:
    transformer_backbone=pickle.load(input_file)
with open(dirGlimps + "transformerA.pickle", "rb") as input_file:
    transformerA=pickle.load(input_file)
with open(dirGlimps + "transformerC.pickle", "rb") as input_file:
    transformerC=pickle.load(input_file)

mypdb = mdt.load(referencetop)

sel = mypdb.topology.select('not type H')
mynoHpdb = mypdb.atom_slice(sel)

mynoHpdb.save(outputtop)

sel = mypdb.topology.select('name =~ "C1."')
myc1pdb = mypdb.atom_slice(sel)

#cg_test_backbone = mdt.load(inputtraj, top=myc1pdb.top)[::10]
cg_test_backbone = mdt.load(inputtraj, top=myc1pdb.top)
dx_test_backbone = dx.ComplementaryDuplex(cg_test_backbone)
lendx=len(dx_test_backbone)

cg_test_backbone

dx_test_backbone.sequence#,dx_test_backbone.sequence==Sequences[simnum-1]

#cg_test_backbone.save("AtomAAAA_skip.mdcrd")

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

fg_predtest.save(outputtraj)

print("Backmapping complete")

# ## Final prediction evaluation

# from mdtraj.geometry.alignment import rmsd_qcp

# sel = fg_predtest.topology.select('name =~ "C1."')
# cg_predtest = fg_predtest.atom_slice(sel)

# 
# rmsd_cg_test=np.zeros(len(cg_test_backbone))
# for i in range(len(cg_test_backbone)):
#     rmsd_cg_test[i]=rmsd_qcp(cg_test_backbone.xyz[i],cg_predtest.xyz[i])

# rmsd_cg_test.mean()*10 #nm -> 10A
