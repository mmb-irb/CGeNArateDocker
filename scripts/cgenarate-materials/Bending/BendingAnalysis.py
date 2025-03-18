path = 'NAFlex/CURVES/'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

twist = pd.read_csv(path + 'NAFlex_canalOut_twist.ser', header=None, sep=r'\s+', index_col=0).T.dropna()
tilt = pd.read_csv(path + 'NAFlex_canalOut_tilt.ser', header=None, sep=r'\s+', index_col=0).T.dropna()
roll = pd.read_csv(path + 'NAFlex_canalOut_roll.ser', header=None, sep=r'\s+', index_col=0).T.dropna()

cum_twist = (twist.cumsum(axis=0))* 2*np.pi /360

rcos = roll * np.cos(cum_twist)
rsin = roll * np.sin(cum_twist)

tcos = tilt * np.cos(cum_twist)
tsin = tilt * np.sin(cum_twist)

bx = (rcos + tsin)
by = (rsin + tcos)

Indices = bx.index

N = len(cum_twist.index)

def bx_nk(n,k):
    if n+k-1 > N:
        return None 
    return bx.loc[n:n+k-1].sum(axis=0)

def by_nk(n,k):
    if n+k-1 > N:
        return None 
    return by.loc[n:n+k-1].sum(axis=0)

def Btot_nk(n,k):
    return np.sqrt(bx_nk(n,k)**2 + by_nk(n,k)**2) 

dflist = []
for n in Indices[:-5]:
    dflist.append((Btot_nk(n,5)))
Btot_n_5 = pd.concat(dflist)

dflist = []
for n in Indices[:-10]:
    dflist.append((Btot_nk(n,10)))
Btot_n_10 = pd.concat(dflist)


dflist = []
for n in Indices[:-5]:
    dflist.append((bx_nk(n,5)))
bx_n_5 = pd.concat(dflist)

dflist = []
for n in Indices[:-10]:
    dflist.append((bx_nk(n,10)))
bx_n_10 = pd.concat(dflist)


dflist = []
for n in Indices[:-5]:
    dflist.append((by_nk(n,5)))
by_n_5 = pd.concat(dflist)

dflist = []
for n in Indices[:-10]:
    dflist.append((by_nk(n,10)))
by_n_10 = pd.concat(dflist)


Btot_1_N = Btot_nk(1,N)
bx_1_N = bx_nk(1,N)
by_1_N = by_nk(1,N)

Btot_n_5.x, Btot_n_5.y = sns.kdeplot(Btot_n_5).lines[-1].get_data()
Btot_n_10.x, Btot_n_10.y = sns.kdeplot(Btot_n_10).lines[-1].get_data()

bx_n_5.x, bx_n_5.y = sns.kdeplot(bx_n_5).lines[-1].get_data()
bx_n_10.x, bx_n_10.y = sns.kdeplot(bx_n_10).lines[-1].get_data()
by_n_5.x, by_n_5.y = sns.kdeplot(by_n_5).lines[-1].get_data()
by_n_10.x, by_n_10.y = sns.kdeplot(by_n_10).lines[-1].get_data()

Btot_1_N.x, Btot_1_N.y = sns.kdeplot(Btot_1_N).lines[-1].get_data()
bx_1_N.x, bx_1_N.y = sns.kdeplot(bx_1_N).lines[-1].get_data()
by_1_N.x, by_1_N.y = sns.kdeplot(by_1_N).lines[-1].get_data()

plt.legend(['Btot_n_5','Btot_n_10','bx_n_5','bx_n_10','by_n_5','by_n_10','Btot_1_N','bx_1_N','by_1_N'])

ind = ['"V{}"'.format(i) for i in range(1, len(Btot_n_5.x)+1)]

file = '"",' + ','.join(ind) + \
    '\n"1",' + ','.join(Btot_n_5.x.astype(str)) + \
    '\n"2",' + ','.join(Btot_n_5.y.astype(str)) + \
    '\n"3",' + ','.join(Btot_n_10.x.astype(str)) + \
    '\n"4",' + ','.join(Btot_n_10.y.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Bending_distribution_total_bending_ensemble.csv', 'w') as f:
    f.write(file)
print(file)

file = '"",' + ','.join(ind) + \
    '\n"1",' + ','.join(bx_n_5.x.astype(str)) + \
    '\n"2",' + ','.join(bx_n_5.y.astype(str)) + \
    '\n"3",' + ','.join(bx_n_10.x.astype(str)) + \
    '\n"4",' + ','.join(bx_n_10.y.astype(str)) + \
    '\n"5",' + ','.join(by_n_5.x.astype(str)) + \
    '\n"6",' + ','.join(by_n_5.y.astype(str)) + \
    '\n"7",' + ','.join(by_n_10.x.astype(str)) + \
    '\n"8",' + ','.join(by_n_10.y.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Bending_distribution_xz_yz_ensemble.csv', 'w') as f:
    f.write(file)
print(file)

file = '"",' + ','.join(ind) + \
    '\n"1",' + ','.join(bx_1_N.x.astype(str)) + \
    '\n"2",' + ','.join(bx_1_N.y.astype(str)) + \
    '\n"3",' + ','.join(by_1_N.x.astype(str)) + \
    '\n"4",' + ','.join(by_1_N.y.astype(str)) + \
    '\n"5",' + ','.join(Btot_1_N.x.astype(str)) + \
    '\n"6",' + ','.join(Btot_1_N.y.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Bending_distribution_total_whole_fiber_ensemble.csv', 'w') as f:
    f.write(file)
print(file)

x_5 = Indices[5:]
x_10 = Indices[10:]

dflist = []
for n in Indices[:-5]:
    dflist.append((bx_nk(n,5)))
bx_a_5 = pd.concat(dflist, axis = 1).mean()

dflist = []
for n in Indices[:-10]:
    dflist.append((bx_nk(n,10)))
bx_a_10 = pd.concat(dflist, axis = 1).mean()


dflist = []
for n in Indices[:-5]:
    dflist.append((by_nk(n,5)))
by_a_5 = pd.concat(dflist, axis = 1).mean()

dflist = []
for n in Indices[:-10]:
    dflist.append((by_nk(n,10)))
by_a_10 = pd.concat(dflist, axis = 1).mean()


dflist = []
for n in Indices[:-5]:
    dflist.append((Btot_nk(n,5)))
Btot_a_5 = pd.concat(dflist, axis = 1).mean()

dflist = []
for n in Indices[:-10]:
    dflist.append((Btot_nk(n,10)))
Btot_a_10 = pd.concat(dflist, axis = 1).mean()


plt.plot(x_5, bx_a_5)
plt.plot(x_10, bx_a_10)
plt.plot(x_5, by_a_5)
plt.plot(x_10, by_a_10)

plt.plot(x_5, Btot_a_5)
plt.plot(x_10, Btot_a_10)

plt.legend(['bx_a_5','bx_a_10','by_a_5','by_a_10','Btot_a_5','Btot_a_10'])

# Plot for 5 values
plt.figure(figsize=(10, 5))
plt.plot(x_5, bx_a_5, label='bx_a_5')
plt.plot(x_5, by_a_5, label='by_a_5')
plt.plot(x_5, Btot_a_5, label='Btot_a_5')
plt.legend()
plt.title('Plot for 5 values')
plt.show()

# Plot for 10 values
plt.figure(figsize=(10, 5))
plt.plot(x_10, bx_a_10, label='bx_a_10')
plt.plot(x_10, by_a_10, label='by_a_10')
plt.plot(x_10, Btot_a_10, label='Btot_a_10')
plt.legend()
plt.title('Plot for 10 values')
plt.show()

ind = ['"V{}"'.format(i) for i in range(1, len(bx_a_5)+1)]

file = '"",' + ','.join(ind) + \
    '\n"1",' + ','.join(x_5.astype(str)) + \
    '\n"2",' + ','.join(bx_a_5.astype(str)) + \
    '\n"3",' + ','.join(x_5.astype(str)) + \
    '\n"4",' + ','.join(by_a_5.astype(str)) + \
    '\n"5",' + ','.join(x_5.astype(str)) + \
    '\n"6",' + ','.join(Btot_a_5.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Individual_bending_xz_yz_ensemble1.csv', 'w') as f:
    f.write(file)
print(file)

file = '"",' + ','.join(ind) + \
    '\n"1",' + ','.join(x_10.astype(str)) + \
    '\n"2",' + ','.join(bx_a_10.astype(str)) + \
    '\n"3",' + ','.join(x_10.astype(str)) + \
    '\n"4",' + ','.join(by_a_10.astype(str)) + \
    '\n"5",' + ','.join(x_10.astype(str)) + \
    '\n"6",' + ','.join(Btot_a_10.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Individual_bending_xz_yz_ensemble2.csv', 'w') as f:
    f.write(file)
print(file)

xz = abs(bx_1_N / Btot_1_N)
yz = abs(by_1_N / Btot_1_N)

xz.x, xz.y = sns.kdeplot(xz).lines[-1].get_data()
yz.x, yz.y = sns.kdeplot(yz).lines[-1].get_data()

plt.legend(["xz","yz"])

ind = ['"V{}"'.format(i) for i in range(1, len(xz.x)+1)]

file = '"",' + ','.join(ind) + \
    '\n"1",' + ','.join(xz.x.astype(str)) + \
    '\n"2",' + ','.join(xz.y.astype(str)) + \
    '\n"3",' + ','.join(yz.x.astype(str)) + \
    '\n"4",' + ','.join(yz.y.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Bending_distribution_perc_total_whole_fiber_ensemble.csv', 'w') as f:
    f.write(file)
print(file)

x_traj = bx.columns

# Unnecessary plot

plt.plot(x_traj, xz)
plt.plot(x_traj, yz)
plt.show()
plt.plot(x_traj, bx_1_N)
plt.plot(x_traj, by_1_N)
plt.plot(x_traj, Btot_1_N)

plt.legend(['xz','yz','Btot_1_N','bx_1_N','by_1_N'])

ind = ['"V{}"'.format(i) for i in range(1, len(x_traj)+1)]

file = '"",' + ','.join(ind) + \
    '\n"",' + ','.join(x_traj.astype(str)) + \
    '\n"d1",' + ','.join(xz.astype(str)) + \
    '\n"",' + ','.join(x_traj.astype(str)) + \
    '\n"d2",' + ','.join(yz.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Bending_total_perc_whole_fiber_alongtraj_ensemble.csv', 'w') as f:
    f.write(file)
print(file)

file = '"",' + ','.join(ind) + \
    '\n"",' + ','.join(x_traj.astype(str)) + \
    '\n"d1",' + ','.join(bx_1_N.astype(str)) + \
    '\n"",' + ','.join(x_traj.astype(str)) + \
    '\n"d2",' + ','.join(by_1_N.astype(str)) + \
    '\n"",' + ','.join(x_traj.astype(str)) + \
    '\n"d3",' + ','.join(Btot_1_N.astype(str))
    #  ','.join(Writhe.apply("{:.03f}".format))
with open('Bending/Bending_total_whole_fiber_alongtraj_ensemble.csv', 'w') as f:
    f.write(file)
print(file)



