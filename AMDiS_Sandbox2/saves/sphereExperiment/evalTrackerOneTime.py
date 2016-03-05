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


ti = 1.7

phi0 = ph1[1][ph1[0]==ti]
phi1 = ph2[1][ph2[0]==ti]
phi2 = ph3[1][ph3[0]==ti]
phi3 = ph4[1][ph4[0]==ti]

if (phi0 < 0): phi0 += 2*pi
if (phi1 < 0): phi1 += 2*pi
if (phi2 < 0): phi2 += 2*pi
if (phi3 < 0): phi3 += 2*pi



#os = 1.0
#hy = 0.2
#yT = arange(os,os+4*hy,hy)
yT=[1.]
p0 = ax.plot(phi0, yT, 'o')
p1 = ax.plot(phi1, yT, 'o')
p2 = ax.plot(phi2, yT, 'o')
p3 = ax.plot(phi3, yT, 'o')
#print plt.getp(p1)

aw= 0.1
arr0 = plt.Arrow(phi0, 0.0, 0, 1, width=aw, color=p0[0].get_color())
arr1 = plt.Arrow(phi1, 0.0, 0, 1, width=aw, color=p1[0].get_color())
arr2 = plt.Arrow(phi2, 0.0, 0, 1, width=aw, color=p2[0].get_color())
arr3 = plt.Arrow(phi3, 0.0, 0, 1, width=aw, color=p3[0].get_color())
ax.add_artist(arr0)
ax.add_artist(arr1)
ax.add_artist(arr2)
ax.add_artist(arr3)

angres = 100
angYs = linspace(yT[0]/10., yT[0]/4. ,4)
lw = 2
plot(linspace(0, phi0, angres), [angYs[0]]*angres, color=p0[0].get_color(), linewidth=lw)
plot(linspace(0, phi1, angres), [angYs[1]]*angres, color=p1[0].get_color(), linewidth=lw)
plot(linspace(0, phi2, angres), [angYs[2]]*angres, color=p2[0].get_color(), linewidth=lw)
plot(linspace(0, phi3, angres), [angYs[3]]*angres, color=p3[0].get_color(), linewidth=lw)



yL=[r'$t_2$']
plt.yticks(yT, yL)

xT=[0, pi/2, pi, 3*pi/2]
xL=[r'$0$', r'$\frac{\pi}{2}$', r'$\pm\pi$', r'$-\frac{\pi}{2}$']
plt.xticks(xT, xL)

xlabel(r'$\varphi$')

ylim([0,yT[0]*1.1])

#circle = plt.Circle((0, 0), os, transform=ax.transData._b, color="gray", alpha=0.5)
#ax.add_artist(circle)

grid(True)
#legend()
show()
