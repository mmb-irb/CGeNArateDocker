# Basic data. DNA residue names (AMBER) and derived lists and dictionaries.
resnamemap = {
              'A': ['DA5', 'DA', 'DA3'],
              'C': ['DC5', 'DC', 'DC3'],
              'G': ['DG5', 'DG', 'DG3'],
              'T': ['DT5', 'DT', 'DT3']
             }
resnames = [v for l in resnamemap.values() for v in l]
rescode = {}
for k, v in resnamemap.items():
    for r in v:
        rescode[r] = k
