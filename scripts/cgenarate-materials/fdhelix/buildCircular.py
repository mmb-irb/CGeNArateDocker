#!/usr/bin/env python
# coding: utf-8

import numpy as np
import sys

if len(sys.argv) != 4: 
	print("Usage: script.py <inputfile> <outputfile> <delta>") 
	sys.exit(1) 

inputfile = sys.argv[1] 
outputfile = sys.argv[2] 
delta_link = sys.argv[3]

#inputfile = 'ATAT.pdb'
#outputfile = 'ATAT2.pdb'
#delta_link=0

canonical_turns=None
target_turns=None

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


# Precomputing angle scaling

xlist=[]
ylist=[]
zlist=[]

rise = 3.4

with open(inputfile, 'r')  as fin:
    for index, line in enumerate(fin):
        if "C1'" in line:
            xlist.append(float(line[46-16:54-16]))
            ylist.append(float(line[46-8:54-8]))
            zlist.append(float(line[46:55]))
print(len(xlist))

zdist=max(zlist)-min(zlist) + rise

x=np.array(xlist)
y=np.array(ylist)
z=np.array(zlist)

n=len(x)
bp=n//2
period = 10.4
if not canonical_turns:
    canonical_turns=round(bp/period)
if not target_turns:
    target_turns=round(canonical_turns+float(delta_link))

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

print(canonical_turns, bp/10.4, target_turns)

turns=(np.sum(allangles[:]))/(2*np.pi)
turnsconnecting=turns+connectangle/(2*np.pi)

angleproportion=(target_turns-1/period)/turns

print(turns, turnsconnecting, (target_turns-1/period)-turns)

xlist=[]
ylist=[]
zlist=[]

with open(inputfile, 'r')  as fin:
    for index, line in enumerate(fin):
        if "ATOM" in line:
            xlist.append(float(line[46-16:54-16]))
            ylist.append(float(line[46-8:54-8]))
            zlist.append(float(line[46:55]))
        else:
            print(line)
print(len(xlist))

# zdist=max(zlist)-min(zlist)

x=np.array(xlist)
y=np.array(ylist)
z=np.array(zlist)

n=len(x)
bp=n//2

origin=originx,originy=0,0

points=np.array((x,y)).T

pointsini=np.array((1,0))
pointsfin=points[[0,-1]]
pointsori=points[[0,-1]]*0+origin

angle_0,angle_last=angles(pointsini,pointsori,pointsfin)

pointsini=points[:-1]
pointsfin=points[1:]
pointsori=points[:-1]*0+origin

allangles=angles(pointsini,pointsori,pointsfin)

alldists = np.sqrt(x**2 + y**2)

# print(z[:100],z[1::30] / zdist,(target_turns-turns))

# newangles = allangles.cumsum() * angleproportion + angle_0
# print(newangles[::30])
newangles = allangles.cumsum() + (z[1:] - min(z)) / zdist * ((target_turns-1/period)-turns)*(2*np.pi) + angle_0
# print(newangles[::30])
#newangles = np.concatenate((newangles + angle_0, newangles[::-1]+angle_last))

newx, newy = x.copy(), y.copy()

newx[1:], newy[1:] = alldists[1:] * np.cos(
    newangles), alldists[1:] * np.sin(newangles)

x, y = newx, newy

r=y+zdist/(2*np.pi)
theta= z/zdist*(2*np.pi)
r,theta

y2, z2 = r*np.cos(theta), r*np.sin(theta)
y2, z2 = y2+zdist/(2*np.pi), z2+zdist/(2*np.pi)
y2, z2

index=-1
with open(inputfile, 'r')  as fin:
    with open(outputfile, 'w')  as fout:
        for i, line in enumerate(fin):
            if "ATOM" in line:
                index+=1
                #print(line.strip())
                line2 = line[:4] + "  {:5d}".format(index+1) + line[11:30]
                line2 += "{:8.2f}".format(x[index]) + "{:8.2f}".format(y2[index]) + "{:8.2f}".format(z2[index]) + line[54:]
                
    #             print(line2.strip())
                fout.writelines(line2)
            else:
                fout.writelines(line)


