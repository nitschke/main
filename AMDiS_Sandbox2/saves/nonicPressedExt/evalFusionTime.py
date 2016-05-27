#!/usr/bin/python

import csv
from pylab import *

fn = "fusionTime.csv"

lineStyles = ['-','--', '-.', ':']
lw = 3;

with open(fn, 'rb') as f:
    reader = csv.DictReader(f, skipinitialspace = True)
    stretch=[];
    press=[];
    ftime=[];
    for row in reader:
        stretch.append(row['stretch'])
        press.append(row['press'])
        ftime.append(row['fusionTime'])
stretch = array(stretch,dtype=float)
press = array(press,dtype=float)
ftime = array(ftime,dtype=float)

fig = plt.figure()
ax1 = fig.add_subplot(111)

#semilogy(stretch, ftime, label='Fusion Time', linestyle="-", marker="o", linewidth=lw)
plot(stretch, ftime, label='Fusion Time', linestyle="-", marker="o", linewidth=lw)

cFusion = 0.635 #(0.83325+0.8335)/2. # lin. interpol.
axvspan(0.0,cFusion, facecolor='0.5', alpha=0.5)
text(0.4, 5, "Not Stable 4 Defects", horizontalalignment='center')
text(1.4, 5, "Stable 4 Defects", horizontalalignment='center')

#xlim(0.0,stretch[-1])
xlim(0.0,2.0)
xlabel('Stretch Factor C')
ylabel('Time t')

locator_params(nbins=10) #don't work with semilogy plot
grid(True)
legend()

# C = (20/7)*B
ax2 = ax1.twiny()
ax2.set_xlabel('Press Factor B')
#ax2.set_xlim(0.0,press[-1])
ax2.set_xlim(0.0,0.7)
#grid(True)

show()
