import numpy as np
import mdtraj as mdt
from . import data

def load(trajfiles, topfile=None, **kwargs):
    """
    Load one or more trajectory files, only keeping DNA atoms.

    """
    if topfile is None:
        if isinstance(trajfiles, list):
            topf = trajfiles[0]
        else:
            topf = trajfiles
    else:
        topf = topfile
    top = mdt.load_topology(topf)
    indices = nucleic(top)
    if 'stride' in kwargs:
        stride = kwargs['stride']
    else:
        stride = 1
    if topfile is None:
        traj = mdt.load(trajfiles, atom_indices=indices, stride=stride)
    else:
        traj = mdt.load(trajfiles, top=topfile, atom_indices=indices, stride=stride)
    return process(traj)

def nucleic(top):
    """
    Returns indices of atoms in traj that correspond to nucleic acid residues

    Arg:
        top: MDTraj topology

    Returns:
        np.array
    """
    n_names = data.resnames
    indices = [a.index for a in top.atoms 
               if a.residue.name in n_names]
    return np.array(indices, dtype=np.int32)

def process(traj):
    """
    Returns a new trajectory for just DNA residues
    
    One chain for each strand

    Args:
        traj: MDTraj trajectory

    Returns:
        MDTraj trajectory
    """
    indices = nucleic(traj.topology)
    nuctraj = traj.atom_slice(np.array(indices))
    # fix cases where multiple strands in same chain:
    if traj.topology.n_chains == 1: # typical if topology comes from prmtop file
        mols = nuctraj.topology.find_molecules()
        cid = np.zeros(nuctraj.n_atoms, dtype=np.int32)
        ic = -1
        for mol in mols:
            ic += 1
            for a in mol:
                cid[a.index] = ic
        atoms, bonds = nuctraj.topology.to_dataframe()
        atoms['chainID'] = cid
        newtop = mdt.Topology().from_dataframe(atoms, bonds)
                
        nuctraj.top = newtop
    return nuctraj

def is_dna(traj):
    """
    Returns True if all residues are DNA residues

    Args:
        traj: MDTraj trajectory

    Returns:
        Bool
    """
    n_names  = data.resnames
    is_dna = True
    for r in traj.topology.residues:
        is_dna = is_dna and r.name in n_names
    return is_dna

def sequences(traj):
    """
    The base sequence in t, one string per strand

    Args:
        traj: MDTraj trajectory

    Returns:
        list of strings
    """
    code = data.rescode
    if not is_dna(traj):
        raise ValueError('Error: trajectory is not pure DNA')
    seqs = []
    for c in traj.topology.chains:
        seq = ''.join([code[r.name] for r in c.residues])
        seqs.append(seq)
    return seqs
