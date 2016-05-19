#!/usr/bin/python

import csv
from pylab import *

fn = "finalEnergy.csv"

lineStyles = ['-','--', '-.', ':']
lw = 3;

with open(fn, 'rb') as f:
    reader = csv.DictReader(f, skipinitialspace = True)
    stretch=[];
    press=[];
    e2D=[];
    e4D=[];
    for row in reader:
        stretch.append(row['stretch'])
        press.append(row['press'])
        e2D.append(row['EnergyOf2Defects'])
        e4D.append(row['EnergyOf4Defects'])
stretch = array(stretch,dtype=float)
press = array(press,dtype=float)
e2D = array(e2D,dtype=float)
e4D = array(e4D,dtype=float)

fig = plt.figure()
ax1 = fig.add_subplot(111)

plot(stretch, e2D, label='2 Defects', linewidth=lw)

cFusion = 0.625 #(0.83325+0.8335)/2. # lin. interpol.
dataFilter = stretch > cFusion
plot(stretch[dataFilter], e4D[dataFilter], '*-',label='4 Defects', linewidth=lw)

axvspan(0.0,cFusion, facecolor='0.5', alpha=0.5)
text(0.4, 16.2, "Not Stable 4 Defects", horizontalalignment='center')
text(1.4, 16.2, "Stable 4 Defects", horizontalalignment='center')

xlim(0.0,2.0)
xlabel('Stretch Factor C')
ylabel('Energy E')

locator_params(nbins=10)
grid(True)
legend()

# C = (20/7)*B
ax2 = ax1.twiny()
ax2.set_xlabel('Press Factor B')
ax2.set_xlim(0.0,0.7)
#grid(True)

show()
