import pandas as pd

df = pd.read_csv('NAFlex/CURVES/NAFlex_canalOut_twist.ser', header=None, sep=r'\s+', index_col=0).T.dropna()

Helical = (df.sum()+df.mean())/360
Lk = round(Helical[1])
print('Helical:', Helical[1], 'Lk:', Lk)
Writhe = (Lk - Helical)

ind = ['"V{}"'.format(i) for i in range(1, len(Writhe)+1)]
integers = ['{}'.format(i) for i in range(1, len(Writhe)+1)]

fileW = '"",' + ','.join(ind) + '\n"1",' + ','.join(integers) + '\n"2",' + ','.join(Writhe.apply("{:.03f}".format))
with open('Circular/Writhe.csv', 'w') as f:
    f.write(fileW)

fileT = '"",' + ','.join(ind) + '\n"1",' + ','.join(integers) + '\n"2",' + ','.join(Helical.apply("{:.03f}".format))
with open('Circular/twist.csv', 'w') as f:
    f.write(fileT)
