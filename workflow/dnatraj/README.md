# DNAtraj

DNAtraj adds a layer above [MDTraj](http://mdtraj.org) designed specifically for
DNA simulation data. It doesn't supply any analysis tools itself, but helps with
the process of loading, editing and transforming raw trajectory data into a format
that can easily be interfaced with other Python-based analysis tools.

## Example:

```python
from dnatraj import duplex as dx
```

Load data from a trajectory. Everything that is not DNA is ignored. The .load() method works just the same way as MDTraj's and supports the same wide range of file formats.


```python
dup = dx.load('seq01_rep6_nv_5_frames.nc', top='repartitioned_seq01.prmtop')
```

"dup" is a ComplementaryDuplex. This has some basic attributes:


```python
print(dup)
print(dup.sequence)
print(len(dup)) # This is the length of the duplex, not the length of the trajectory!
```

    <ComplementaryDuplex GCACTAGGCTAGCCTAGTGC with 5 frames>
    GCACTAGGCTAGCCTAGTGC
    20


`ComplementaryDuplex`es can be sliced, when it makes sense:


```python
core = dup[2:18] # remove terminal GC base pairs
print(core)
print(dup[2:-2])
try:
    print(dup[2:18:4])
except:
    print('slices must be contiguous')
try:
    print(dup[14:12])
except:
    print('slicing backwards is not allowed')
```

    <ComplementaryDuplex ACTAGGCTAGCCTAGT with 5 frames>
    <ComplementaryDuplex ACTAGGCTAGCCTAGT with 5 frames>
    slices must be contiguous
    slicing backwards is not allowed


`ComplementaryDuplex`es can be inverted (strands swapped):


```python
acta = core[:4]
print(acta)
print(acta.invert())
```

    <ComplementaryDuplex ACTA with 5 frames>
    <ComplementaryDuplex TAGT with 5 frames>


`ComplementaryDuplex`es have a `.where()` method that extracts all subduplexes that match the given pattern.


```python
for cta in core.where('CTA'):
    print(cta)
```

    <ComplementaryDuplex CTA with 5 frames>
    <ComplementaryDuplex CTA with 5 frames>
    <ComplementaryDuplex CTA with 5 frames>
    <ComplementaryDuplex CTA with 5 frames>
    <ComplementaryDuplex CTA with 5 frames>
    <ComplementaryDuplex CTA with 5 frames>


By default this searches both Watson and Crick strands, but this can be changed:


```python
for cta in core.where('CTA', mode='W'):
    print(cta)
```

    <ComplementaryDuplex CTA with 5 frames>
    <ComplementaryDuplex CTA with 5 frames>
    <ComplementaryDuplex CTA with 5 frames>


`ComplementaryDuplex`es can be saved in a variety of formats. As well as the trajectory, you can save a compatible "topology" file in .pdb, .gro, or AMBER .prmtop format. If you choose .pdb or .gro then the file will contain original atom and residue numbers, but if you don't want this, you can renumber eveything (atoms and residues). Prmtop files created by the program are minimalistic - topology is OK but forcefield parameters are bogus. It should be good enough to be acceptable to, e.g. ccptraj though.


```python
ctas = core.where('CTA')
ctas[0].save('cta_instance_0.nc', top='cta.prmtop')
ctas[4].renumber()
ctas[4].save('cta_instance_4_renumbered.nc')
```

`ComplementaryDuplex`es can be concatenated along the time axis (as long as they have the same sequence):


```python
cta = dx.tstack(ctas)
print(cta)
cta.save('cta_6.nc')
```

    <ComplementaryDuplex CTA with 30 frames>


They can also be concatenated along the sequence axis (as long as they have the same number of frames):


```python
a = dup[2:5]
b = dup[5:7]
c = dx.sstack([a, b])
print(a)
print(b)
print(c)
```

    <ComplementaryDuplex TAG with 5 frames>
    <ComplementaryDuplex GC with 5 frames>
    <ComplementaryDuplex TAGGC with 5 frames>

**Note**: it is the user's responsibility to ensure the coordinates of the resulting duplex are sensible.


When you want to do some analysis, the underlying MDTraj trajectory object is accessible:


```python
print(cta.traj)
```

    <mdtraj.Trajectory with 30 frames, 191 atoms, 6 residues, and unitcells>



## Author:
[Charlie Laughton](mailto:charles.laughton@nottingham.ac.uk)
