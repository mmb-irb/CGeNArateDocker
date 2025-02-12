#!/usr/bin/env python3
#
# Fake an Amber prmtop format file from a PDB file.
#
# Topology info should be legitimate, but all parameters (force constants,
# charges, etc.) will be bogus.
# Should contain just enough info to be acceptable to programs that need
# a prmtop file to support loading trajectory file data,  e.g. Chimera
#
# Works as both a command line tool and importable module.
#
# Command line:
#
#    fakeparm.py <pdbfile> <prmtopfile>
#
#
# As module:
#
#  from fakeparm import top2prmtopfile
#
#  top = traj.topology # an MDTraj topology
#  prmtopfile = 'fake.prmtop'
#  top2prmtopfile(top, prmtopfile)
#  
import mdtraj as mdt
import numpy as np
import networkx as nx
from operator import itemgetter
from argparse import ArgumentParser

def fortwrite(f, itemfmt, itemsperline, data):
    '''
    little routine to do fortran-style output formatting
    '''
    indx = 0
    if len(data) == 0:
        f.write('\n')
        return
    while indx < len(data):
        nitems = min(itemsperline, len(data[indx:]))
        fmt = itemfmt * nitems + '\n'
        f.write(fmt.format(*data[indx:indx+nitems]))
        indx += nitems

def top2prmtopfile(top, prmtopfile):

    # Use networkx methods to find bonds, angles and dihedrals
    g = nx.Graph()
    g.add_nodes_from(top.atoms)
    g.add_edges_from(top.bonds)
    atoms = list(top.atoms)
    bonds = []
    angles = []
    dihedrals = []
    exclusions = []
    nexcl = []
    for a in atoms:
        nex = 0
        i = a.index
        sp = nx.single_source_shortest_path(g, a, 4)
        nex = 0
        for t in sp:
            p = sp[t]
            if p[-1].index > i:
                exclusions.append(t.index + 1)
                nex += 1
                if len(p) == 2:
                    bonds.append(p)
                elif len(p) == 3:
                    angles.append(p)
                elif len(p) == 4:
                    dihedrals.append(p)
        if nex == 0:
            exclusions.append(0)
            nex = 1
        nexcl.append(nex)

    # Process data to create required prmtop parameters
    ITITL = 'Fake prmtop file' 
    NATOM = top.n_atoms
    elements = list(set([a.element.symbol for a in atoms]))
    NTYPES = len(elements)
    bonh = [b for b in bonds if b[0].element.symbol == 'H' or b[1].element.symbol == 'H']
    NBONH = len(bonh)
    bona = [b for b in bonds if b[0].element.symbol != 'H' and b[1].element.symbol != 'H']
    MBONA = len(bona)
    angh = [a for a in angles if a[0].element.symbol == 'H' or a[2].element.symbol == 'H']
    NTHETH = len(angh)
    anga = [a for a in angles if a[0].element.symbol != 'H' and a[2].element.symbol != 'H']
    MTHETA = len(anga)
    dih = [d for d in dihedrals if d[0].element.symbol == 'H' or d[3].element.symbol == 'H']
    NPHIH = len(dih)
    dia = [d for d in dihedrals if d[0].element.symbol != 'H' and d[3].element.symbol != 'H']
    MPHIA = len(dia)
    NHPARM = 0
    NPARM = 0
    NUMEX = nexcl
    INB = exclusions
    NNB = len(INB)
    LBRES = ['{:4s}'.format(r.name) for r in top.residues]
    NRES = len(LBRES)
    NBONA = MBONA
    NTHETA = MTHETA
    NPHIA = MPHIA
    NUMBND = min(1, NBONH + MBONA)
    NUMANG = min(1, NTHETH + MTHETA)
    NPTRA = min(1, NPHIH + MPHIA)
    NATYP = NTYPES
    NPHB = 0
    IFPERT = 0
    NBPER = 0
    NGPER = 0
    NDPER = 0
    MBPER = 0
    MGPER = 0
    MDPER = 0
    IFBOX = 0 # NB: no periodic box info is assumed
    NMXRS = 0
    for r in top.residues:
        NMXRS = max(NMXRS, len(list(r.atoms)))
    IFCAP = 0
    NUMEXTRA = 0
    NCOPY = 0
    IGRAPH = ['{:4s}'.format(a.name) for a in atoms]
    CHARGE = [0.0] * NATOM # All atom charges are zero
    ATNUM = [a.element.atomic_number for a in atoms]
    AMASS = [a.element.mass for a in atoms]
    IAC = [elements.index(a.element.symbol) + 1 for a in atoms]
    ICO = []
    for i in range(NTYPES):
        for j in range(NTYPES):
            ICO.append(i + j + 1)
    IPRES = [r._atoms[0].index + 1 for r in top.residues]
    # Bogus forcefield parameter values follow:
    RK = [0.0] * NUMBND
    REQ = [3.0] * NUMBND
    TK = [0.0] * NUMANG
    TEQ = [0.0] * NUMANG
    PK = [0.0] * NPTRA
    PN = [1] * NPTRA
    PHASE = [0.0] * NPTRA
    SOLTY = [0] * NATYP
    CN1 = [1.0] * (NTYPES * (NTYPES+1) // 2)
    CN2 = CN1
    # Create the bond, angle, and dihedral lists
    IBH = []
    bonh = [[b[0].index, b[1].index] for b in bonh]
    bonh.sort(key=itemgetter(0, 1))
    for b in bonh:
        IBH.append(b[0] * 3)
        IBH.append(b[1] * 3)
        IBH.append(1)
    IB = []
    bona = [[b[0].index, b[1].index] for b in bona]
    bona.sort(key=itemgetter(0, 1))
    for b in bona:
        IB.append(b[0] * 3)
        IB.append(b[1] * 3)
        IB.append(1)
    ITH = []
    angh = [[a[0].index, a[1].index, a[2].index] for a in angh]
    angh.sort(key=itemgetter(1, 0, 2))
    for a in angh:
        ITH.append(a[0] * 3)
        ITH.append(a[1] * 3)
        ITH.append(a[2] * 3)
        ITH.append(1)
    IT = []
    anga = [[a[0].index, a[1].index, a[2].index] for a in anga]
    anga.sort(key=itemgetter(1, 0, 2))
    for a in anga:
        IT.append(a[0] * 3)
        IT.append(a[1] * 3)
        IT.append(a[2] * 3)
        IT.append(1)
    IPH = []
    dih = [[d[0].index, d[1].index, d[2].index, d[3].index] for d in dih]
    dih.sort(key=itemgetter(1, 2, 0, 3))
    for d in dih:
        IPH.append(d[0] * 3)
        IPH.append(d[1] * 3)
        IPH.append(d[2] * 3)
        IPH.append(d[3] * 3)
        IPH.append(1)
    IP = []
    dia = [[d[0].index, d[1].index, d[2].index, d[3].index] for d in dia]
    dia.sort(key=itemgetter(1, 2, 0, 3))
    for d in dia:
        IP.append(d[0] * 3)
        IP.append(d[1] * 3)
        IP.append(d[2] * 3)
        IP.append(d[3] * 3)
        IP.append(1)
    # A few final bits of data:
    ISYMBL = ['{:4s}'.format(a.element.symbol) for a in atoms]
    ITREE = ['BLA '] * NATOM
    JOIN = [0] * NATOM
    IROTAT = [0] * NATOM
    TYPE = 'RADIUS_SET'
    RBORN = [2.0] * NATOM
    FS = [1.0] * NATOM

    # Create the prmtop file:
    with open(prmtopfile, 'w') as f:
        f.write('%VERSION VERSION_STAMP = V0001.000 DATE = 05/22/06  12:10:21\n')
        f.write('%FLAG TITLE\n%FORMAT(20A4)\n')
        f.write('{:80s}\n'.format(ITITL))
        
        f.write('%FLAG POINTERS\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, [NATOM, NTYPES, NBONH, MBONA, NTHETH, MTHETA,
                                   NPHIH, MPHIA, NHPARM, NPARM,
                                   NNB, NRES, NBONA, NTHETA, NPHIA, NUMBND, 
                                   NUMANG, NPTRA, NATYP, NPHB,
                                   IFPERT, NBPER, NGPER, NDPER, MBPER, MGPER, 
                                   MDPER, IFBOX, NMXRS, IFCAP,
                                   NUMEXTRA, NCOPY])

        f.write('%FLAG ATOM_NAME\n%FORMAT(20a4)\n')
        fortwrite(f, '{:4s}', 20, IGRAPH)
        
        f.write('%FLAG CHARGE\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, CHARGE)
        
        #f.write('%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n')
        #fortwrite(f, '{:8d}', 10, ATNUM)
        
        f.write('%FLAG MASS\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, AMASS)
        
        f.write('%FLAG ATOM_TYPE_INDEX\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IAC)
        
        f.write('%FLAG NUMBER_EXCLUDED_ATOMS\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, NUMEX)
        
        f.write('%FLAG NONBONDED_PARM_INDEX\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, ICO)
        
        f.write('%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\n')
        fortwrite(f, '{:4s}', 20, LBRES)
        
        f.write('%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IPRES)
        
        f.write('%FLAG BOND_FORCE_CONSTANT\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, RK)
        
        f.write('%FLAG BOND_EQUIL_VALUE\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, REQ)
        
        f.write('%FLAG ANGLE_FORCE_CONSTANT\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, TK)
        
        f.write('%FLAG ANGLE_EQUIL_VALUE\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, TEQ)
        
        f.write('%FLAG DIHEDRAL_FORCE_CONSTANT\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, PK)
        
        f.write('%FLAG DIHEDRAL_PERIODICITY\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, PN)
        
        f.write('%FLAG DIHEDRAL_PHASE\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, PHASE)
        
        f.write('%FLAG SOLTY\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, SOLTY)
        
        f.write('%FLAG LENNARD_JONES_ACOEF\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, CN1)
        
        f.write('%FLAG LENNARD_JONES_BCOEF\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, CN2)
        
        f.write('%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IBH)
        
        f.write('%FLAG BONDS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IB)
        
        f.write('%FLAG ANGLES_INC_HYDROGEN\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, ITH)
        
        f.write('%FLAG ANGLES_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IT)
        
        f.write('%FLAG DIHEDRALS_INC_HYDROGEN\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IPH)
        
        f.write('%FLAG DIHEDRALS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IP)
        
        f.write('%FLAG EXCLUDED_ATOMS_LIST\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, INB)

        f.write('%FLAG HBOND_ACOEF\n%FORMAT(5E16.8)\n\n')
        f.write('%FLAG HBOND_BCOEF\n%FORMAT(5E16.8)\n\n')
        f.write('%FLAG HBCUT\n%FORMAT(5E16.8)\n\n')

        f.write('%FLAG AMBER_ATOM_TYPE\n%FORMAT(20a4)\n')
        fortwrite(f, '{:4s}', 20, ISYMBL)
        
        f.write('%FLAG TREE_CHAIN_CLASSIFICATION\n%FORMAT(20a4)\n')
        fortwrite(f, '{:4s}', 20, ITREE)
        
        f.write('%FLAG JOIN_ARRAY\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, JOIN)
        
        f.write('%FLAG IROTAT\n%FORMAT(10I8)\n')
        fortwrite(f, '{:8d}', 10, IROTAT)
        
        f.write('%FLAG RADIUS_SET\n%FORMAT(1A80)\n{:80s}\n'.format(TYPE))
        
        f.write('%FLAG RADII\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, RBORN)
        
        f.write('%FLAG SCREEN\n%FORMAT(5E16.8)\n')
        fortwrite(f, '{:16.8e}', 5, FS)
        

if __name__ == '__main__':
    parser = ArgumentParser(description='Create a fake prmtop file from a pdb file')
    parser.add_argument('pdbfile', help='name of input pdb file')
    parser.add_argument('prmtopfile', help='name of prmtop file to create')
    args = parser.parse_args()
    traj = mdt.load(args.pdbfile)
    top = traj.topology
    top2prmtopfile(top, args.prmtopfile)

