from re import T
import mdtraj as mdt
import numpy as np
import dnatraj
from dnatraj import sequencetools
from dnatraj.fakeparm import top2prmtopfile
from mdtraj import Trajectory as MDTrajectory

class Duplex(object):
    def __init__(self, source):
        if isinstance(source, MDTrajectory):     
            self.traj = dnatraj.process(source)
        elif isinstance(source, (Duplex, ComplementaryDuplex)):
                self.traj = source.traj
        if not _is_duplex(self.traj):
            raise TypeError('Error: trajectory is not a duplex')
        self._sequences = dnatraj.sequences(self.traj)
        self._len = len(self._sequences[0])
        
    @property
    def sequences(self):
        return self._sequences
    
    @property
    def sequence(self):
        return self._sequences
    
    def __len__(self):
        return self._len

    def renumber(self):
        ir = 0
        ia = 0
        for r in self.traj.topology.residues:
            ir += 1
            r.resSeq = ir
            for a in r.atoms:
                ia += 1
                a.serial = ia

    def save(self, filename, top=None):
        if top is not None:
            try:
                self.traj[0].save(top)
            except:
                top2prmtopfile(self.traj.topology, top)
        self.traj.save(filename)

    def copy(self):
        return self.__class__(self.traj[:])
       
    def slice(self, key):
        if isinstance(key, int):
            c = _subduplex(self.traj, key, 1)
        elif isinstance(key, slice):
            if key.step is not None:
                raise ValueError('Error: only contiguous slices are allowed')
            if key.start is None:
                start = 0
            else:
                start = key.start
            if key.stop is None:
                stop = self._len
            else:
                stop = key.stop
            if stop < 0:
                stop = self._len + stop
            if stop < start:
                raise IndexError('Error: negative slice ranges not permitted.')
            c = _subduplex(self.traj, start, stop-start)
        else:
            raise ValueError('Error, key must be int or slice')
        return ComplementaryDuplex(c)
    
    def invert(self):
        inv = _invert(self.traj)
        return self.__class__(inv)
    
    def __getitem__(self, key):
        return self.slice(key)
    
    def __str__(self):
        return '<{} {} with {} frames>'.format(self.__class__.__name__, self.sequence, self.traj.n_frames)
    
    def __repr__(self):
        return '<{} at {}>'.format(self.__class__.__name__, id(self))
    
class ComplementaryDuplex(Duplex):
    def __init__(self, source):
        super().__init__(source)
        if not _is_wc(self.traj):
            raise TypeError('Error: trajectory is not a complementary duplex')
        
    @property
    def sequence(self):
        return self._sequences[0]
    
    def where(self, subseq, mode='WC'):
        subtrajs = _subduplexes(self.traj, subseq, mode=mode)
        subcds = [ComplementaryDuplex(st) for st in subtrajs]
        return subcds

def load(traj, top=None, **kwargs):
    """
    Load a duplex

    Args:
        traj: str, name of trajectory file
        top: str, name of topoplogy file

    Returns:
        Duplex or ComplementaryDuplex
    """
    n = dnatraj.load(traj, topfile=top, **kwargs)
    if not _is_duplex(n):
        raise TypeError('Error: trajectory is not a duplex')
    if _is_wc(n):
        return ComplementaryDuplex(n)
    else:
        return Duplex(n)

def _is_duplex(traj):
    """
    Test if a trajectory is for double-stranded DNA

    There must be two strands of equal length, but they do not have to be
    complementary

    Args:
        traj: MDTraj trajectory

    Returns:
        Bool
    """
    seqs = dnatraj.sequences(traj)
    return len(seqs) == 2 and (len(seqs[0]) == len(seqs[1]))

def _is_wc(traj):
    """
    Test if a trajectory is a Watson-Crick paired duplex.

    Note the test is for complementary sequences only, not
    H-bonding patterns.

    Args:
        traj: MDTraj trajectory

    Returns:
        Bool
    """
    result = _is_duplex(traj)
    if result:
        seqs = dnatraj.sequences(traj)
        result = seqs[0] == sequencetools.complement(seqs[1])
    return result

def _sequence(traj):
    """
    Base sequence of the duplex (Watson strand)

    Args:
        traj: MDTraj trajectory

    Returns:
        str, base sequence of Watson strand
    """
    if not _is_wc(traj):
        raise TypeError('Error: trajectory is not a WC duplex')
    seqs = dnatraj.sequences(traj)
    return seqs[0]

def _sequences(traj):
    """
    Base sequence of the duplex (both strands)

    Args:
        traj: MDTraj trajectory

    Returns:
        (str, str), base sequencase of Watson and Crick strands
    """
    if not _is_duplex(traj):
        raise TypeError('Error: trajectory is not a duplex')
    seqs = dnatraj.sequences(traj)
    return seqs

def _invert(traj):
    """
    Invert the two strands in a DNA duplex

    Args:
        traj: MDTraj trajectory

    Returns:
        MDTraj trajectory
    """
    if not _is_duplex(traj):
        raise TypeError('Error: trajectory is not a duplex')
    s1 = traj.topology.select('chainid 0')
    s2 = traj.topology.select('chainid 1')
    it = traj.topology.subset(s2)
    it = it.join(traj.topology.subset(s1))
    xyz = traj.xyz[:, np.concatenate((s2, s1))]
    return mdt.Trajectory(xyz, it, time=traj.time, 
                          unitcell_lengths=traj.unitcell_lengths, 
                          unitcell_angles=traj.unitcell_angles)

def _subduplex(traj, start, length):
    """
    Extract a subsequence from a duplex

    Args:
        traj: MDTraj trajectory
        start: int, index with offset of 1st base to extract
        length: int, number of base pairs to extract

    Returns:
        MDTraj trajectory
    """
    if not _is_duplex(traj):
        raise TypeError('Error: trajectory is not a duplex')
    seq = _sequence(traj)    
    l = len(seq)
    if start >= l or (start + length) > l:
        raise IndexError('Error: duplex is only {} bases long'.format(l))
    s1 = start
    e1 = start + length - 1
    s2 = 2 * l - e1 - 1
    e2 = s2 + length - 1
    sel = 'resid {} to {} or resid {} to {}'.format(s1, e1, s2, e2)
    indices = traj.topology.select(sel)
    return traj.atom_slice(indices)
    

def _subsequence_locations(traj, subseq):
    """
    Returns the locations of the given subsequence in t

    Args:
        traj: MDTraj trajectory
        subseq: string, the base sequence to search for

    Returns:
        list with two lists of offsets into the associated strand where the
            first base in the subsequence is found.
    """
    if not _is_duplex(traj):
        raise TypeError('Error: trajectory is not a duplex')
    seq = dnatraj.sequences(traj)
    locs = []
    for strand in seq:
        slocs = []
        l = strand.find(subseq)
        while l > -1:
            slocs.append(l)
            l = strand.find(subseq, l+1)
        locs.append(slocs)
    return locs  

def _subduplexes(traj, subseq, mode='WC'):
    """
    Returns a list of trajectories, one for each occurrence of
    the subsequence

    Args:
        traj: MDtraj trajectory
        subseq: str, base sequence to match
        mode: str, whether to search both strands ("WC"), or just the
              first ("W") or second ("C")

    Returns:
        list of MDtraj trajectories
    """
    locs = _subsequence_locations(traj, subseq)
    l = len(subseq)
    sts = []
    if 'W' in mode:
        st0 = [_subduplex(traj, s, l) for s in locs[0]]
        sts += st0
    if 'C' in mode:
        traj = _invert(traj)
        st1 = [_subduplex(traj, s, l) for s in locs[1]]
        sts += st1
    return sts

def _wc_hbonds(traj, heavy=False):
    """
    Return WC Hbond data

    Args:
        traj: MDTraj trajectory
        heavy: Bool, if True, report heavy atom-heavy atom distances.

    Returns:
        dict, keys are H-bond identifiers, values are arrays of distances,
              one per snapshot in traj.
    """
    hbdefs = {
        'A':['resid {} and name N1 or resid {} and name H3',
            'resid {} and name H61 or resid {} and name O4'],
        'T':['resid {} and name O4 or resid {} and name H61',
            'resid {} and name H3 or resid {} and name N1'],
        'C':['resid {} and name N3 or resid {} and name H1',
            'resid {} and name O2 or resid {} and name H21',
            'resid {} and name H41 or resid {} and name O6'],
        'G':['resid {} and name H1 or resid {} and name N3',
            'resid {} and name H21 or resid {} and name O2',
            'resid {} and name O6 or resid {} and name H41']
    }
    
    hbdefsH = {
        'A':['resid {} and name N1 or resid {} and name N3',
            'resid {} and name N6 or resid {} and name O4'],
        'T':['resid {} and name O4 or resid {} and name N6',
            'resid {} and name N3 or resid {} and name N1'],
        'C':['resid {} and name N3 or resid {} and name N1',
            'resid {} and name O2 or resid {} and name N2',
            'resid {} and name N4 or resid {} and name O6'],
        'G':['resid {} and name N1 or resid {} and name N3',
            'resid {} and name N2 or resid {} and name O2',
            'resid {} and name O6 or resid {} and name N4']
    }
    if not _is_wc(traj):
        raise TypeError('Error: trajectory is not a WC duplex')
    seqs = dnatraj.sequences(traj)
    nbp = len(seqs[0])
    hbonds = {}
    for iw in range(nbp):
        ic = 2 * nbp - iw - 1
        if heavy:
            hbps = hbdefsH[traj.topology.residue(iw).name[1]]
        else:
            hbps = hbdefs[traj.topology.residue(iw).name[1]]
        for hbp in hbps:
            sel = hbp.format(iw, ic)
            wc = traj.topology.select(sel)
            key = '{}:{}'.format(traj.topology.atom(wc[0]), traj.topology.atom(wc[1]))
            hbonds[key] = wc
    z = [hbonds[k] for k in hbonds]
    hbdists = mdt.compute_distances(traj, z)
    wchb = {k:v for k, v in zip(hbonds.keys(), hbdists.T)}
    return wchb

def _tstack(trajs):
    """
    Stack a list of trajectories along the time axis
    
    Args:
        trajs: list of MDTraj trajectories
        
    Returns:
    
    MDTraj trajectory
    """

    if not isinstance(trajs, list):
        trajs = [trajs]
    if len(trajs) == 1:
        result = trajs[0]
    else:
        result = trajs[0].join(trajs[1:])
    return result

def tstack(dups):
    """
    Stack a list of Duplexes along the time axis
    
    Args:
        dups: list of Duplexes
        
    Returns:
        Duplex
    """

    if not isinstance(dups, list):
        dups = [dups]
    for dup in dups:
        if not isinstance(dup, (Duplex, ComplementaryDuplex)):
            raise TypeError('Error: arguments must be duplexes')

    trajs = [d.traj for d in dups]
    return Duplex(_tstack(trajs))

def _sstack(trajs):
    """
    Stack a list of duplex trajectories along the sequence axis
    
    Args:
        trajs: list of MD Trajectories
        
    Returns:

    MDTraj trajectory
    
    """

    if not isinstance(trajs, list):
        trajs = [trajs]
    for t in trajs:
        if not _is_duplex(t):
            raise TypeError('Error - trajectory is not a duplex')
    nf = trajs[0].n_frames
    for t in trajs[1:]:
        if t.n_frames != nf:
            raise ValueError('Error - trajectories not all same length')
    
    tstacked = None
    for t in trajs:
        sel = t.topology.select('chainid 0')
        tw = t.atom_slice(sel)
        if tstacked is None:
            tstacked = tw
        else:
            tstacked = tstacked.stack(tw)
    
    for t in trajs[::-1]:
        sel = t.topology.select('chainid 1')
        tc = t.atom_slice(sel)
        tstacked = tstacked.stack(tc)
    df, bonds = tstacked.topology.to_dataframe()
    nr = tstacked.topology.n_residues
    last_a = tstacked.topology.select(f'resid {(nr//2)-1}')[-1]
    df.loc[:last_a+1, 'chainID'] = 0
    df.loc[last_a+1:, 'chainID'] = 1
    top = mdt.Topology().from_dataframe(df, bonds)
    tstacked.topology = top
    
    return tstacked

def sstack(dups):
    """
    Stack a list of duplexes along the sequence axis
    
    Args:
       dups: list of duplexes
       
    Returns:
       Duplex
    """

    if not isinstance(dups, list):
        dups = [dups]
    for dup in dups:
        if not isinstance(dup, (Duplex, ComplementaryDuplex)):
            raise TypeError('Error: arguments must be duplexes')

    trajs = [d.traj for d in dups]
    dup = Duplex(_sstack(trajs))
    try:
        dup = ComplementaryDuplex(dup)
    except TypeError:
        pass
    return dup
