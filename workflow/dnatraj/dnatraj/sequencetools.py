# tools that work on strings, not trajectories
def is_dna(seq):
    """
    Checks sequence is DNA
    """
    bases = set([b for b in seq])
    return bases <= set(('A', 'C', 'G', 'T'))

def complement(seq):
    """
    Returns complement sequence to seq
    """
    if not is_dna(seq):
        raise ValueError('Error: sequence contains letter other than ACGT')
    t = seq.translate(seq.maketrans('ACGT','TGCA'))
    return t[::-1]
 
def ry(seq):
    """
    Converts base sequence to RY sequence
    """
    if not is_dna(seq):
        raise ValueError('Error: sequence contains letter other than ACGT')
    t = seq.translate(seq.maketrans('ACGT','RYRY'))
    return t

