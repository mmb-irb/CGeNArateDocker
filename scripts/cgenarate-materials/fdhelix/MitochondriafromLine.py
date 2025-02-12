#!/usr/bin/env python
# coding: utf-8

# In[1]:


inputfile = 'ATAT.pdb'
outputfile = 'ATAT2.pdb'


# In[2]:


canonical_turns=None
delta_link=0
target_turns=None


# In[3]:


import numpy as np


# In[4]:


def angles(a, b, c):

    ba = a - b
    bc = c - b
    ac = c - a

    ba /= np.expand_dims(np.linalg.norm(ba, axis=-1), axis=-1)
    bc /= np.expand_dims(np.linalg.norm(bc, axis=-1), axis=-1)

    cosine_angle = (ba * bc).sum(axis=-1)
    #     cosine_angle=np.minimum(cosine_angle,1)
    angle = np.arccos(cosine_angle)
    
    rotate=np.array(((0,-1),(1,0)))

    sign = np.sign((bc @ rotate.T * ac).sum(axis=-1))
    return angle*sign


# In[5]:


xlist=[]
with open(inputfile, 'r')  as fin:
    for index, line in enumerate(fin):
        xlist.append(float(line[46-16:54-16]))
        
ylist=[]
with open(inputfile, 'r')  as fin:
    for index, line in enumerate(fin):
        ylist.append(float(line[46-8:54-8]))
        
zlist=[]
with open(inputfile, 'r')  as fin:
    for index, line in enumerate(fin):
        zlist.append(float(line[46:55]))


# In[6]:


zdist=max(zlist)-min(zlist)

x=np.array(xlist)
y=np.array(ylist)
z=np.array(zlist)

n=len(x)
bp=n//2


# In[7]:


if not canonical_turns:
    canonical_turns=round(bp/10.4)
if not target_turns:
    target_turns=canonical_turns+delta_link


# In[8]:


origin=originx,originy=0,0

points=np.array((x,y)).T

pointsini=np.array((1,0))
pointsfin=points[[0,-1]]
pointsori=points[[0,-1]]*0+origin

angle_0,angle_last=angles(pointsini,pointsori,pointsfin)

pointsini=points[:-bp-1]
pointsfin=points[1:bp]
pointsori=points[:-bp-1]*0+origin

allangles=angles(pointsini,pointsori,pointsfin)

pointsini=points[-bp-1:-bp]
pointsfin=points[:1]
pointsori=points[:1]*0+origin

connectangle=angles(pointsini,pointsori,pointsfin)[0]


# In[9]:


turns=(np.sum(allangles[:]))/(2*np.pi)
turnsconnecting=turns+connectangle/(2*np.pi)

angleproportion=target_turns/turnsconnecting

turns, turnsconnecting


# In[10]:


alldists = np.sqrt(x**2 + y**2)

newangles = allangles.cumsum() * angleproportion 
newangles = np.concatenate((newangles + angle_0, newangles[::-1]+angle_last))

newx, newy = x.copy(), y.copy()

newx[1:-1], newy[1:-1] = alldists[1:-1] * np.cos(
    newangles), alldists[1:-1] * np.sin(newangles)


# In[11]:


x, y = newx, newy


# In[12]:


r=y+zdist/(2*np.pi)
theta= z/zdist*(2*np.pi)
r,theta


# In[13]:


y2, z2 = r*np.cos(theta), r*np.sin(theta)
y2, z2 = y2+zdist/(2*np.pi), z2+zdist/(2*np.pi)
y2, z2


# In[14]:


with open(inputfile, 'r')  as fin:
    with open(outputfile, 'w')  as fout:
        for index, line in enumerate(fin):
#             print(line.strip())
            line2 = line[:4] + "  {:5d}".format(index+1) + line[11:30]
            line2 += "{:8.2f}".format(x[index]) + "{:8.2f}".format(y2[index]) + "{:8.2f}".format(z2[index]) + line[54:]
            
#             print(line2.strip())
            fout.writelines(line2)

