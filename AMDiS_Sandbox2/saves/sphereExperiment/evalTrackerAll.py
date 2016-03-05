#!/usr/bin/python

import csv
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fn = 'sphereExtremeValues.csv'

with open(fn, 'rb') as f:
    reader = csv.DictReader(f, skipinitialspace = True)
    t = [ndarray((0,1),dtype=double)]
    phi = [ndarray((0,1),dtype=double)]
    for row in reader:
        t = append(t, double(row["t"]))
        x = double(row["x"])
        y = double(row["y"])
        val = x / sqrt(x*x + y*y)
        if y >= 0.0 :
            phi = append(phi, arccos(val))
        else :
            phi = append(phi, -pi + arccos(-val))

offFilter = t>0.001
t = t[offFilter]        
phi = phi[offFilter]
ph1 = [[t[0],phi[0]]]
ph2 = [[t[1],phi[1]]]
ph3 = [[t[2],phi[2]]]
ph4 = [[t[3],phi[3]]]
eps = 0.1
for i in arange(4,len(t)):
    if (abs(ph1[-1][1] - phi[i]) < eps):
        ph1.append([t[i], phi[i]])
    if (abs(ph2[-1][1] - phi[i]) < eps):
        ph2.append([t[i], phi[i]])
    if (abs(ph3[-1][1] - phi[i]) < eps):
        ph3.append([t[i], phi[i]])
    if (abs(ph4[-1][1] - phi[i]) < eps):
        ph4.append([t[i], phi[i]])

ph1 = transpose(array(ph1));        
ph2 = transpose(array(ph2));        
ph3 = transpose(array(ph3));        
ph4 = transpose(array(ph4));        

        
#scatter(t, phi)
rcParams['lines.linewidth'] = 4
#ls = ['-', '--', ':', '-.']
ls = ['-', '-', '-', '-']

fig = plt.figure()
ax = fig.gca(polar=True)


t0 = 0.002
t1 = 0.5
t2 = 1.7
t3 = 4.5


os= 5.0
yT = [t0+os, t1+os, t2+os, t3+os]
#ax.scatter(phi,t+os)
ax.plot(ph1[1],ph1[0]+os, linestyle=ls[0])
ax.plot(ph2[1],ph2[0]+os, linestyle=ls[1])
ax.plot(ph3[1],ph3[0]+os, linestyle=ls[2])
ax.plot(ph4[1],ph4[0]+os, linestyle=ls[3])

yL=[r'$t_0$', r'$t_1$', r'$t_2$', r'$t_3$']
plt.yticks(yT, yL)


xT=[0, pi/2, pi, 3*pi/2]
xL=[r'$0$', r'$\frac{\pi}{2}$', r'$\pm\pi$', r'$-\frac{\pi}{2}$']
plt.xticks(xT, xL)

xlabel(r'$\varphi$')
ylim(0,t[-1]+os)

circle = plt.Circle((0, 0), os, transform=ax.transData._b, color="lightgray", alpha=1)
ax.add_artist(circle)

grid(True)
#legend()
show()
