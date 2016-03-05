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
        
        
#scatter(t, phi)

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

fig = plt.figure()
ax = fig.gca(polar=True)


t0 = 0.002
t1 = 0.5
t2 = 1.7
t3 = 4.5

phi0 = concatenate((ph1[1][ph1[0]==t0], ph1[1][ph1[0]==t1], ph1[1][ph1[0]==t2], ph1[1][ph1[0]==t3]))
phi1 = concatenate((ph2[1][ph2[0]==t0], ph2[1][ph2[0]==t1], ph2[1][ph2[0]==t2], ph2[1][ph2[0]==t3]))
phi2 = concatenate((ph3[1][ph3[0]==t0], ph3[1][ph3[0]==t1], ph3[1][ph3[0]==t2], ph3[1][ph3[0]==t3]))
phi3 = concatenate((ph4[1][ph4[0]==t0], ph4[1][ph4[0]==t1], ph4[1][ph4[0]==t2], ph4[1][ph4[0]==t3]))


os = 1.0
hy = 0.2
yT = arange(os,os+4*hy,hy)
ax.plot(phi0, yT[arange(len(phi0))],'o')
ax.plot(phi1, yT[arange(len(phi1))],'o')
ax.plot(phi2, yT[arange(len(phi2))],'o')
ax.plot(phi3, yT[arange(len(phi3))],'o')


yL=[r'$t_0$', r'$t_1$', r'$t_2$', r'$t_3$']
plt.yticks(yT, yL)

xT=[0, pi/2, pi, 3*pi/2]
xL=[r'$0$', r'$\frac{\pi}{2}$', r'$\pm\pi$', r'$-\frac{\pi}{2}$']
plt.xticks(xT, xL)

xlabel(r'$\varphi$')
ylim(0,yT[-1]+hy/5)

circle = plt.Circle((0, 0), os, transform=ax.transData._b, color="gray", alpha=0.5)
ax.add_artist(circle)

grid(True)
#legend()
show()
