#!/usr/bin/env python
# coding: utf-8

# We need to know
# 
# 1. DNA sequence
# 1. TF position
# 1. TF length
# 1. TF pdb

import os
import numpy as np
from scipy.spatial.distance import cdist
import sys
import argparse

#string = '''3 input/1vfc.pdb 12
#30 input/1wtr.pdb 7
#60 input/3zhm.pdb 8
#80 input/1bnz.pdb 7'''

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

# Create the input folder if it doesn't exist
if not os.path.exists('input'):
    os.makedirs('input')

# Arguments
#Sequence_FILE='input/Proteins.seq'
Sequence_FILE='inputSequence.txt'
data_outfile="input/Proteins.txt"
# pdb_list = ['input/' + pdb for pdb in str.split()[1::3]]
pdb_list = string.split()[1::3]
print(pdb_list)
dyads = np.array(string.split()[::3],int)
dyads_len = np.array(string.split()[2::3],int)

# Auxiliary functions
def exyz_to_q(ex, ey, ez):
    """ taken from LAMMPS source code, converts ex,ey,ez frame vectors to quaternion orientations """

    # squares of quaternion components

    q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0)
    q1sq = q0sq - 0.5 * (ey[1] + ez[2])
    q2sq = q0sq - 0.5 * (ex[0] + ez[2])
    q3sq = q0sq - 0.5 * (ex[0] + ey[1])

    q = np.array([0.0, 0.0, 0.0, 0.0])
    # some component must be greater than 1/4 since they sum to 1
    # compute other components from it

    if q0sq >= 0.25:
        q[0] = np.sqrt(q0sq)
        q[1] = (ey[2] - ez[1]) / (4.0 * q[0])
        q[2] = (ez[0] - ex[2]) / (4.0 * q[0])
        q[3] = (ex[1] - ey[0]) / (4.0 * q[0])
    elif q1sq >= 0.25:
        q[1] = np.sqrt(q1sq)
        q[0] = (ey[2] - ez[1]) / (4.0 * q[1])
        q[2] = (ey[0] + ex[1]) / (4.0 * q[1])
        q[3] = (ex[2] + ez[0]) / (4.0 * q[1])
    elif q2sq >= 0.25:
        q[2] = np.sqrt(q2sq)
        q[0] = (ez[0] - ex[2]) / (4.0 * q[2])
        q[1] = (ey[0] + ex[1]) / (4.0 * q[2])
        q[3] = (ez[1] + ey[2]) / (4.0 * q[2])
    elif q3sq >= 0.25:
        q[3] = np.sqrt(q3sq)
        q[0] = (ex[1] - ey[0]) / (4.0 * q[3])
        q[1] = (ez[0] + ex[2]) / (4.0 * q[3])
        q[2] = (ez[1] + ey[2]) / (4.0 * q[3])

    norm = np.linalg.norm(q)
    q = q / norm

    return q

def q_to_exyz(q):
    ex = (q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3],
          2.0 * (q[1] * q[2] + q[0] * q[3]),
          2.0 * (q[1] * q[3] - q[0] * q[2]))

    ey = (2.0 * (q[1] * q[2] - q[0] * q[3]),
          q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3],
          2.0 * (q[2] * q[3] + q[0] * q[1]))

    ez = (2.0 * (q[1] * q[3] + q[0] * q[2]),
          2.0 * (q[2] * q[3] - q[0] * q[1]),
          q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3])

    return np.array(ex), np.array(ey), np.array(ez)

def get_com(data):
    return np.mean(data, axis=0)

def get_diam(data):
    return np.max(cdist(data, data, 'euclidean'))

def get_nucl_dna(DNA_FILE):
    
    data=list(open(DNA_FILE))
    
    coords=[]
    for line in data:
        if (" C1' " in line and "ATOM" in line):
            sline=line.split()
            x=np.array(sline[6:9],dtype=float)
            coords.append(x)
            # print(line)

    return np.array(coords)

def get_core(CA_FILE):
    data=list(open(CA_FILE))
    
    coords=[]
    chains=[]
    for line in data:
        if (" CA " in line and "ATOM" in line):
            # print(line)
            sline=line.split()
            x=np.array(sline[6:9],dtype=float)
            chain=sline[4]
            chains.append(chain)
            coords.append(x)

    return np.array(coords),chains

def fit_core(data):
    # first compute the reference frame vector then convert to quat

    # origin is the com
    com = get_com(data)

    # # x axis points to specific beads
    # x1=113
    # x2=600
    # x=(data[x1,:]+data[x2,:])*0.5
    # xv = x-com
    
    # #zef 
    # z1=466
    # z2=953

    # zef1=data[z1,:]
    # zef2=data[z2,:]

    # zefv=zef1-zef2

    # yv=-np.cross(xv,zefv)

    # zv=np.cross(xv,yv)

    # xv=unit_vec(xv)
    # yv=unit_vec(yv)
    # zv=unit_vec(zv)

    # q=exyz_to_q(xv,yv,zv)

    # We don't care about orientation.
    q = exyz_to_q((1,0,0),(0,1,0),(0,0,1))

    return com,q

def find_naked_regions(N, dyads_position, dyads_len):
    linker_begins = np.append([0],dyads_position+dyads_len)
    linker_ends = np.append(dyads_position-2,[N//2-1])

    linker_regions = [(b,e) for b,e in zip(linker_begins, linker_ends)]
    
    return linker_regions

def linker_from_scratch(linker_N):
    W=np.full((linker_N,3),[2.308, -5.333, -0.452])
    C=np.full((linker_N,3),[2.308,  5.333,  0.452])

    RISE=3.4
    TWIST=35*np.pi/180.0

    W[:,2]+=np.arange(linker_N)*RISE
    C[:,2]+=np.arange(linker_N)*RISE

    TWIST_array=np.arange(linker_N)*TWIST

    rot=np.array([[np.cos(TWIST_array),-np.sin(TWIST_array)],[np.sin(TWIST_array),np.cos(TWIST_array)]])
    rot=np.transpose(rot, (2,1,0))

    W[:,:2]=np.einsum('ijk,ij->ik', rot,W[:,:2])
    C[:,:2]=np.einsum('ijk,ij->ik', rot,C[:,:2])

    dnalinker=np.concatenate((W,C[::-1]),axis=0)
    
    return dnalinker

def find_isometry(target_coords, starting_coords):
    
    if target_coords.shape[0]<4 and target_coords.shape[1]!=3 and         starting_coords.shape[0]<4 and starting_coords.shape[1]!=3:
        print('Not a viable isometry')
        sys.exit(1)
    target2=target_coords-target_coords.mean(axis=0)
    starting2=starting_coords-starting_coords.mean(axis=0)

    H=starting2.T@target2
    U,S,V=np.linalg.svd(H)
    R=V.T@U.T
    
    Trans=-R@starting_coords.mean(axis=0)+target_coords.mean(axis=0)
    
    return R,Trans

def apply_isometry(isometry, starting_coords):
    R,Trans = isometry
    new_coords=(R@starting_coords.T).T+Trans
    return new_coords

def find_dna_isometry(current_dna, next_dna):
    n1=len(current_dna)
    current_tetrad=current_dna[n1//2-2:n1//2+2]
    next_tetrad=next_dna[[0,1,-2,-1]]

    # Create glue
    glue=linker_from_scratch(4)
    glue_start=glue[[0,1,-2,-1]]
    iso1=find_isometry(current_tetrad,glue_start)
    glue_t=apply_isometry(iso1, glue)

    # Transform next dna to nucl
    glue_t_end=glue_t[8//2-2:8//2+2]
    iso2=find_isometry(glue_t_end,next_tetrad)

    return iso2

def apply_dna_isometry(isometry, starting_coords):
    return apply_isometry(isometry, starting_coords)

def join_dna(dna_list):
    new_dna=dna_list[0]
    for next_dna in dna_list[1:]:
        new_dna=join_2_dna(new_dna,next_dna)
    return new_dna

def join_2_dna(current_dna, next_dna):
    new_dna = np.concatenate((current_dna[:len(current_dna)//2], next_dna ,current_dna[len(current_dna)//2:]))
    return new_dna

def number_tetramers(sequence):
    
    def seqtrans(string):
        mydict={'A':'T','C':'G','G':'C','T':'A'}
        mymap='ACGT'.maketrans(mydict)
        return string[::-1].translate(mymap)

    NTetramers={}

    unique=[a+b+c+d for a in 'ACGT' for b in 'ACGT' for c in 'ACGT' for d in 'ACGT']

    for i, element in enumerate(unique):
        NTetramers[element]=i

    Tetramers=[None]*len(NTetramers)
    i=1


    for element in NTetramers:
        if Tetramers[NTetramers[seqtrans(element)]] is None: # and Tetramers[NTetramers[element]] is None
            Tetramers[NTetramers[element]]=i

            i=i+1
        else:
            Tetramers[NTetramers[element]]=-Tetramers[NTetramers[seqtrans(element)]]
            
    sequence='G'+sequence[:len(sequence)//2]+'C'
    
    SegmentedNSequence=[sequence[i: i+4] for i in range(len(sequence)-3)]
    SegmentedSequence=np.vectorize(lambda x: Tetramers[NTetramers[x]])(np.array(SegmentedNSequence))
    shift=(SegmentedSequence>0)*2-1
    SegmentedSequence=SegmentedSequence*shift
    
    tetramers = [(i,j) for i,j in zip(SegmentedSequence,shift)]
                 
    return tetramers

def make_data(cores, all_dna,total_ids,total_types,total_molids, DC, DD, data_outfile,sequence):
    
    atoms=[]
    ells=[]
    bonds=[]
    angles=[]
    
    total_xs=[]

    
    natomtypes=5#len(np.unique(total_types))

    # first cores
    DCX=DC[0]
    DCY=DC[1]
    DCZ=DC[2]
    
    i=0
    for core in cores:
        com = core[0]
        total_xs.append(com)
        q = core[1]
        atoms.append(str(total_ids[i])+ " "+str(total_types[i])+" " + str(com[0])+" "+str(com[1])+" "+str(com[2])+" 1 1 "+str(total_molids[i])+" 0\n")
        ells.append(str(total_ids[i]) + " "+DCX+" "+DCY+" "+DCZ+" "+str(q[0])+" "+str(q[1])+" "+str(q[2])+" "+str(q[3])+"\n")
        i+=1
        
        
# WHAT ARE THE 1 1 n 0
# atom-ID atom-type x y z ellipsoidflag density molecule-ID q
    # dna
    tetramers=number_tetramers(sequence)
    parselist=[ 2,  1,  3,  5,  4,  6,  8,  7,  9, 11, 10, 12]
    parselist2=[ 4,  3,  2,  1]
    
    
    q=np.array((1,0,0,0))
    dids=[]
    for dna in all_dna:
        com=dna
        total_xs.append(com)
        DDX=DDY=DDZ=str(DD[total_types[i]])
        atoms.append(str(total_ids[i])+ " "+str(total_types[i])+" " + str(com[0])+" "+str(com[1])+" "+str(com[2])+" 1 1 "+str(total_molids[i])+" 0\n")
        ells.append(str(total_ids[i])+ " "+DDX+" "+DDY+" "+DDZ+" "+str(q[0])+" "+str(q[1])+" "+str(q[2])+" "+str(q[3])+"\n")
        dids.append(total_ids[i])
        i+=1
    
    # dna-dna bonds
    
    b=1
    n=len(dids)
    for i in range(1,len(dids)//2):
        tetnum,parse=tetramers[i-1]
        parse=(parse==-1)
#         print('Index', i)
        
        auxlist=[None]*4
        for j, k in zip((i-2,i-1,i,i+1),range(4)):
            if j>=0 and j<n//2:
                auxlist[k]=dids[j]
#             else:
#                 print('Skipped {}'.format(j)) 
        id_1,id0,id1,id2=auxlist
                
        auxlist=[None]*4
        for j, k in zip((i-2,i-1,i,i+1),range(4)):
            if j>=0 and j<n//2:
                auxlist[k]=dids[n-j-1]
#                 print('Skipped n-{j}-1')
        idn_1,idn0,idn1,idn2=auxlist
        
        bond = [b,17,id_1,id2]        
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1   
        bond = [b,17,idn_1,idn2]
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1 
#         bond = [b,17,id_1,id1]
#         if None not in bond:
#             bonds.append(" ".join(str(x) for x in bond)+"\n")
#             b+=1   
#         bond = [b,17,idn_1,idn1]
#         if None not in bond:
#             bonds.append(" ".join(str(x) for x in bond)+"\n")
#             b+=1 
#         /*Bond, Bond-*/
#         BondForceij(F, r, DTetramers[TNums[i]][0], KTetramers[TNums[i]][0], i, i + 1);
#         BondForceij(F, r, DTetramers[TNums[i]][1], KTetramers[TNums[i]][1], n - 1 - i, n - 2 - i);
        bond = [b,1,id0,id1]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        bond = [b,2,idn0,idn1]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1

#         /*3m, 2mA, 2mB, 1m*/
#         BondForceij(F, r, DTetramers[TNums[i]][6], KTetramers[TNums[i]][6], i - 1, n - i - 3);
#         BondForceij(F, r, DTetramers[TNums[i]][7], KTetramers[TNums[i]][7], i - 1, n - i - 2);
#         BondForceij(F, r, DTetramers[TNums[i]][8], KTetramers[TNums[i]][8], i, n - i - 3);
#         BondForceij(F, r, DTetramers[TNums[i]][9], KTetramers[TNums[i]][9], i, n - i - 2);
        bond = [b,3,id_1,idn2]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        bond = [b,4,id_1,idn1]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        bond = [b,5,id0,idn2]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        bond = [b,6,id0,idn1]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        
#         /*bpA, bpB*/
#         BondForceij(F, r, DTetramers[TNums[i]][10], KTetramers[TNums[i]][10], i, n - i - 1);
#         BondForceij(F, r, DTetramers[TNums[i]][11], KTetramers[TNums[i]][11], i + 1, n - i - 2);
        
        
        if None not in bond:
            bond = [b,7,id0,idn0]
            bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        if i==1 and None not in bond:
            bond = [b,7,id0,idn0]
            bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
            
        if None not in bond:
            bond = [b,8,id1,idn1]
            bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        if i==len(dids)//2-1 and None not in bond:
            bond = [b,8,id1,idn1]
            bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1

#         /*1p, 2pA, 2pB, 3p*/
#         BondForceij(F, r, DTetramers[TNums[i]][12], KTetramers[TNums[i]][12], i + 1, n - i - 1);
#         BondForceij(F, r, DTetramers[TNums[i]][13], KTetramers[TNums[i]][13], i + 1, n - i);
#         BondForceij(F, r, DTetramers[TNums[i]][14], KTetramers[TNums[i]][14], i + 2, n - i - 1);
#         BondForceij(F, r, DTetramers[TNums[i]][15], KTetramers[TNums[i]][15], i + 2, n - i);
        bond = [b,9,id1,idn0]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        bond = [b,10,id1,idn_1]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        bond = [b,11,id2,idn0]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        bond = [b,12,id2,idn_1]
        bond[1] = tetnum * 100+ bond[1]*(1-parse)+parselist[bond[1]-1]*parse
        if None not in bond:
            bonds.append(" ".join(str(x) for x in bond)+"\n")
            b+=1
        
    
#     /*LONG*/
#     for (i = 0; i < n / 2 - 5; i++) {
#         /*5m, 4m*/
#         BondForceij(F, r, DLong[0], KLong[0], i, n - i - 6);
#         BondForceij(F, r, DLong[1], KLong[1], i, n - i - 5);
#     }
#     /*4m*/
#     i = n / 2 - 5;
#     BondForceij(F, r, DLong[1], KLong[1], i, n - i - 5);
    
#     /*4p*/
#     i = 4;
#     BondForceij(F, r, DLong[2], KLong[2], i, n - i + 3);
#     for (i = 5; i < n / 2; i++) {
#         /*4p, 5p*/
#         BondForceij(F, r, DLong[2], KLong[2], i, n - i + 3);
#         BondForceij(F, r, DLong[3], KLong[3], i, n - i + 4);
#     }
    for i in range(1,len(dids)//2-3):
         
        id0=dids[i-1]
        id4=dids[i+3]
        idn0=dids[n-(i-1)-1]
        idn4=dids[n-(i+3)-1]
        
        bond = [b,13,id0,idn4]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1
        bond = [b,14,id4,idn0]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1
        
        bond = [b,17,id0,id4]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1   
        bond = [b,17,idn0,idn4]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1 
        
    for i in range(1,len(dids)//2-4):
         
        id0=dids[i-1]
        id5=dids[i+4]
        idn0=dids[n-(i-1)-1]
        idn5=dids[n-(i+4)-1]
        
        bond = [b,15,id0,idn5]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1
        bond = [b,16,id5,idn0]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1
        
        bond = [b,17,id0,id5]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1
        bond = [b,17,idn0,idn5]
        bonds.append(" ".join(str(x) for x in bond)+"\n")
        b+=1
        
    nbondtypes=136*100+12
    
    thing=[bond.split()[0] for bond in bonds]
    print('1 is good:',np.unique(thing,return_counts=True)[1].max())
    
    # dna-dna angles
    a=1
    n=len(dids)
    for i in range(1,len(dids)//2):
        tetnum,parse=tetramers[i-1]
        parse=(parse==-1)
        
        auxlist=[None]*4
        for j, k in zip((i-2,i-1,i,i+1),range(4)):
            if j>=0 and j<n//2:
                auxlist[k]=dids[j]
#             else:
#                 print('Skipped {}'.format(j)) 
        id_1,id0,id1,id2=auxlist
                
        auxlist=[None]*4
        for j, k in zip((i-2,i-1,i,i+1),range(4)):
            if j>=0 and j<n//2:
                auxlist[k]=dids[n-j-1]
#                 print('Skipped n-{j}-1')
        idn_1,idn0,idn1,idn2=auxlist

    #         /*AngleA, AngleB, Angle-A, Angle-B*/
#         AngleForceijk(F, r, DTetramers[TNums[i]][2], KTetramers[TNums[i]][2], i - 1, i, i + 1);
#         AngleForceijk(F, r, DTetramers[TNums[i]][3], KTetramers[TNums[i]][3], i, i + 1, i + 2);
#         AngleForceijk(F, r, DTetramers[TNums[i]][4], KTetramers[TNums[i]][4], n - i, n - 1 - i, n - 2 - i);
#         AngleForceijk(F, r, DTetramers[TNums[i]][5], KTetramers[TNums[i]][5], n - 1 - i, n - 2 - i, n - 3 - i);
        angle = [a,1,id_1,id0,id1]
        angle[1] = tetnum * 10+ angle[1]*(1-parse)+parselist2[angle[1]-1]*parse
        if None not in angle:
            angles.append(" ".join(str(x) for x in angle)+"\n")
            a+=1
        angle = [a,2,id0,id1,id2]
        angle[1] = tetnum * 10+ angle[1]*(1-parse)+parselist2[angle[1]-1]*parse
        if None not in angle:
            angles.append(" ".join(str(x) for x in angle)+"\n")
            a+=1
        angle = [a,3,idn_1,idn0,idn1]
        angle[1] = tetnum * 10+ angle[1]*(1-parse)+parselist2[angle[1]-1]*parse
        if None not in angle:
            angles.append(" ".join(str(x) for x in angle)+"\n")
            a+=1
        angle = [a,4,idn0,idn1,idn2]
        angle[1] = tetnum * 10+ angle[1]*(1-parse)+parselist2[angle[1]-1]*parse
        if None not in angle:
            angles.append(" ".join(str(x) for x in angle)+"\n")
            a+=1
        
    thing=[bond.split()[0] for bond in angles]
    print('1 is good:',np.unique(thing,return_counts=True)[1].max())
    
    nangletypes=136*10+4


    total_xs=np.array(total_xs)
    xmin=np.min(total_xs[:,0])-1500
    ymin=np.min(total_xs[:,1])-1500
    zmin=np.min(total_xs[:,2])-500
    
    xmax=np.max(total_xs[:,0])+1500
    ymax=np.max(total_xs[:,1])+1500
    zmax=np.max(total_xs[:,2])+500
    
    data_file = open(data_outfile, "w")
    data_file.write("#LAMMPS data file\n")
    data_file.write("\n")
    data_file.write(str(len(atoms)) + " atoms\n")
    data_file.write(str(len(ells)) + " ellipsoids\n")
    data_file.write(str(len(bonds)) + " bonds\n")
    data_file.write(str(len(angles)) + " angles\n")
    data_file.write("0 dihedrals\n")
    data_file.write("0 impropers\n")
    data_file.write("\n")
    data_file.write(str(natomtypes) +" atom types\n")
    data_file.write(str(nbondtypes) +" bond types\n")
    data_file.write(str(nangletypes) +" angle types\n")
    data_file.write("\n")
    data_file.write(str(xmin)+" "+str(xmax)+" xlo xhi\n")
    data_file.write(str(ymin)+" "+str(ymax)+" ylo yhi\n")
    data_file.write(str(zmin)+" "+str(zmax)+" zlo zhi\n")
    data_file.write("\n")
    data_file.write("Atoms\n")
    data_file.write("\n")
    for line in atoms:
        data_file.write(line)
    data_file.write("\n")
    data_file.write("Ellipsoids\n")
    data_file.write("\n")
    for line in ells:
        data_file.write(line)
    data_file.write("\n")
    data_file.write("Bonds\n")
    data_file.write("\n")
    for line in bonds:
        data_file.write(line)
    data_file.write("\n")
    data_file.write("Angles\n")
    data_file.write("\n")
    for line in angles:
        data_file.write(line)
    data_file.write("\n")
    
    data_file.close()

def make_CGeNArate_Nucleosome(Sequence_FILE, dyads, dyads_len, Nucleosome_FILE=None,                               Different_Nucleosome_FILES=None,                               data_outfile=None, single =True):
    
    DC=["56", "56", "40"]
    
    DD={'A':12.,'C':11.,'G':11.9,'T':10.9}
  
    def seqtrans(string):
        mydict={'A':'T','C':'G','G':'C','T':'A'}
        mymap='ACGT'.maketrans(mydict)
        return string[::-1].translate(mymap)
    
    #We open the sequence File, if it has additional lines we exit
    sequence=list(open(Sequence_FILE))
    if len(sequence) > 1:
        print('There are extra lines in sequence file')
        sys.exit(1)
    sequence=sequence[0].split()[0]
    print(len(sequence))

    if single:
        sequence += seqtrans(sequence)
    print(sequence)
    
    #We check for even number of residues
    N=len(sequence)
    if N%2==1:
        print('Odd number of residues, something is missing')
        sys.exit(1)
        
    #We observe whether the sequence is complementary.
    #If there are > 10% mismatches, there was probably a mistake, and we exit
    #If there are fewer, we inform the user of the mismatches and continue

    W_seq=sequence[:N//2]
    C_seq=sequence[N//2:]
    if seqtrans(C_seq)!=W_seq:
        print('The sequence is not complementary')    
        
        mismatches=np.where([ W_seq[i] != seqtrans(C_seq)[i] for i in range(len(W_seq)) ])[0]
        
        if len(mismatches) > len(W_seq)//10:
            print('There is a {}% of mismatches, probably not intended behaviour'.
                  format(len(mismatches)/len(W_seq)*100))
            sys.exit(1)
            
        for mismatch in mismatches:
            print('Mismatch in Base-pair {}, local context:\n {}\n {}'.
                  format(mismatch,W_seq[max(mismatch-3,0):min(mismatch+4,N//2)],
                         C_seq[N//2-1-max(mismatch-3,0):N//2-1-min(mismatch+4,N//2):-1]))    
        print('Make sure to provide a parametrisation for the relevant Tetramers with mismatches')
        
    ###    
        
    # We get the numeric positions of the dyads from the input
    if type(dyads) == str:
        dyads_position = dyads.split(',')
    else:
        dyads_position = dyads
    
    if dyads_position[0]=='':
        dyads_position=np.array([], dtype=int)
    else:
        dyads_position = np.array(dyads_position, dtype=int)
    
    # If they are too close we exit, also checking they fit within the sequence
    # dyads_diff=np.diff(dyads_position)
    
    # to_exit=False
    # if len(dyads_diff)>0:
    #     if dyads_diff.min() < 147:
    #         print('Some dyads are too close:')
    #         for i in np.where(dyads_diff<147)[0]:
    #             print('Nucleosome {} and {}, with dyads at {} and {}, are at distance {} < 147'.
    #                   format(i,i+1,dyads_position[i],dyads_position[i+1],dyads_diff[i]))
    #         to_exit=True
            
    # if len(dyads_position)>0:
    #     if (dyads_position[-1] > N//2-74) or (dyads_position[0] < 0+74):
    #         print('Dyads should be between 74 and {}'.format(N//2 - 74))
    #         to_exit=True
    # if to_exit:
    #     sys.exit(1)
    print('Assuming TF are far away')
    ###
    
    #If Nucleosome is constant, we distribute it accross all dyads
    if Different_Nucleosome_FILES == None:
        Different_Nucleosome_FILES = [Nucleosome_FILE]*len(dyads_position)
        print("This does not happen, different files")
    
    #If Nucleosome is not constant, we check for same amount of dyads, if not we exit
    elif len(Different_Nucleosome_FILES) != len(dyads_position):
        print('Nucleosome Files {} and dyads {} don\'t match'.format(
            len(Different_Nucleosome_FILES),len(dyads_position)))
        print(Different_Nucleosome_FILES,dyads_position)
        sys.exit(1)

        
    ###
    
    nucl_dna = []
    core_com = []
    core_q = []
    diams = []
    #Read Nucleosomes one by one, creating the corresponding structure
    
    for i, Nucleosome_FILE in enumerate(Different_Nucleosome_FILES):
        # get data from AA nucleosome
        ref_nucl_dna=get_nucl_dna(Nucleosome_FILE) # list of x,q pairs

        core,chains=get_core(Nucleosome_FILE)
        print(core.shape,chains)
        # replace core atoms with the minimal bead
        ref_core_com,ref_core_q = fit_core(core)
        
        diam = get_diam(core)

        nucl_dna.append(ref_nucl_dna)
        core_com.append(ref_core_com)
        core_q.append(ref_core_q)
        diams.append(diam)        
        
    
    #adapt *trim* to create 'de novo' linker DNA, and concatenate with nucleosome


    # program stages:
    # 1. Find which linker dna sequences we need to create, and create it

    # 2. replicate nucleosomes to desired size

    # 3. remap to 5bp resolution using simple method


    # 1.
    print("Creating desired linker DNA")
    
    dna_to_create = find_naked_regions(N, dyads_position,dyads_len)
    linker_dna = [linker_from_scratch(dna_i[1]-dna_i[0]+1) 
                  for dna_i in dna_to_create]
    
    print('Linkers in the following index ranges:', dna_to_create)
#     print([a.shape[0]/2 for a in linker_dna])
#     print(all_nucl_dna,all_core_com,all_core_q)
    assert(len(linker_dna) == len(dyads_position)+1)


    #2.
    print("Placing all DNA Particles consecutively")
    
    all_dna=[linker_dna[0]]
    all_cores=[]
    all_dna_molids=[[0]*len(linker_dna[0])]

    for i in range(len(dyads_position)):

        current_linker_dna = all_dna[-1]
        next_nucl_dna = nucl_dna[i]
        next_linker_dna = linker_dna[i+1]
        next_core_com = core_com[i]
        next_core_q = core_q[i]
        
        # Add nucleosome DNA
        isometry = find_dna_isometry(current_linker_dna, next_nucl_dna)
        next_nucl_dna_t = apply_dna_isometry(isometry, next_nucl_dna)
        all_dna.append(next_nucl_dna_t)
        all_dna_molids.append([i+1]*len(next_nucl_dna_t))
        
        # Add linker DNA
        current_linker_dna = all_dna[-1]
        
        isometry2 = find_dna_isometry(current_linker_dna, next_linker_dna)
        next_linker_dna_t = apply_dna_isometry(isometry2, next_linker_dna)
        all_dna.append(next_linker_dna_t)
        all_dna_molids.append([0]*len(next_linker_dna_t))
        
        # transform a core #FAKE
        T=np.zeros((4,4))
        T[:3,:3],T[:3,3]=isometry

        M=np.zeros((4,4))
        r=next_core_com
        q=next_core_q
        ex,ey,ez=q_to_exyz(q)
        M[:3,0]=ex
        M[:3,1]=ey
        M[:3,2]=ez
        M[:3,3]=r
        M[3,3]=1

        new_M = np.matmul(T,M)
        new_r=new_M[:3,3]

        newex = new_M[:3,0]
        newey = new_M[:3,1]
        newez = new_M[:3,2]
        new_q = exyz_to_q(newex,newey,newez)

#         return T, isometry, M, new_M, new_r, next_nucl_dna

        all_cores.append((new_r,new_q))
        
    print("Joining all DNA Particles")
    dna=join_dna(all_dna)
    # print(len(dna),len(sequence))
    assert(len(dna) == len(sequence))
    
    # correctly assign atom ttypes and ids
    ids = [i for i in range(1,len(dna)+len(all_cores)+1)]

    my_type_dict={'A':2,'C':3,'G':4,'T':5}
    dna_types=[my_type_dict[res] for res in sequence]
    types = [1 for i in range(len(all_cores))] + dna_types
    
    dna_molids = list(join_dna(all_dna_molids))
    molids = [i for i in range(1, len(all_cores)+1)] + dna_molids

    print(len(ids),len(types),len(molids))
#     print(all_cores,dna.shape)

    # 5. write the output file and make bond topology
    DD={my_type_dict[key]: value for key,value in DD.items()}
    # make_data(all_cores,dna,ids,types,molids, DC, DD, data_outfile,sequence)
    
    # print("Wrote "+data_outfile)
    return all_cores,dna,ids,types,molids, DC, DD, data_outfile, sequence, diams

#TETRAMER file something weird first execution?

a=make_CGeNArate_Nucleosome(Sequence_FILE,dyads,dyads_len,data_outfile=data_outfile, Different_Nucleosome_FILES=pdb_list)

text= "all_cores,dna,ids,types,molids, DC, DD, data_outfile, sequence, diams".split(',')
all_cores,dna,ids,types,molids, DC, DD, data_outfile, sequence, diams = a
molids_chain = (np.array(molids)+65).view('U2')

if ids[-1]>10000:
	raise ValueError("Too many residues")

lines = [None] * (len(types)+1)
lines[0] = "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1"

for i,core in enumerate(all_cores,1 + len(dna)):
	lines[i] = "ATOM   {:4d}  NC  NN  {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           C".format(
		i,molids_chain[i-1 - len(dna)],i,*core[0])

for i,(residue,letter) in enumerate(zip(dna,sequence),1):
	lines[i] = "ATOM   {:4d}  C1' D{}  {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           C".format(
		i,letter,molids_chain[i-1 + len(all_cores)],i,*residue)

with open(data_outfile[:-3]+'pdb', 'w') as f:
	f.write('\n'.join(lines) + '\n')

def next_letter(letter):
    return chr(ord(letter) + 1)

lines = [None] * (len(types)+1)
lines[0] = "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1"

for i,core in enumerate(all_cores,1 + len(dna)):
	lines[i] = "ATOM   {:4d}  NC  NN  {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           C".format(
		i,next_letter(molids_chain[i-1 - len(dna)]),i,*core[0])

for i,(residue,letter) in enumerate(zip(dna,sequence),1):
    if i <= len(sequence)//2:
        lines[i] = "ATOM   {:4d}  C1' D{}  {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           C".format(
            i,letter,'A',i,*residue)
    else:
        lines[i] = "ATOM   {:4d}  C1' D{}  {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           C".format(
            i,letter,'B',i,*residue)

with open(data_outfile[:-4]+'Chains.pdb', 'w') as f:
	f.write('\n'.join(lines) + '\n')

lines = [str(len(all_cores)),]
for diam in (diams):
    lines.append(str(diam))

with open('data.in', 'w') as f:
	f.write('\n'.join(lines) + '\n')

