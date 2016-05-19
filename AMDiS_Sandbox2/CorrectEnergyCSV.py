import csv
import sys
import os
from pylab import *
import matplotlib.pyplot as plt

#fn = "output/sphereEnergies.csv"
fn = sys.argv[1]

lineStyles = ['-','--', '-.', ':']
lw = 3;

name = ["Time", "Div", "Rot", "B2", "Norm", "Full"]
n = len(name)
last = 1000000000
with open(fn, 'rb') as f:
    reader = csv.DictReader(f, skipinitialspace = True)
    x = n*[ndarray((0,1),dtype=double)]
    k = 0
    for row in reader:
      for i in arange(n):
        x[i] = append(x[i], double(row[name[i]]))
      if k == last: break
      k = k+1

x[3] = x[3]/2.
x[5] -= x[3] 

nn = len(x[0])
with open(os.path.basename(fn),'w') as f:
    f.write(name[0])
    for i in arange(1,n):
        f.write(',' + name[i])
    for k in arange(0,nn):
        f.write('\n' + str(x[0][k]))
        for i in arange(1,n):
            f.write(',' + str(x[i][k]))
    

#for i in arange(1,n):
#  plot(x[0], x[i], label=name[i], linewidth=lw)
#
#    
#
#grid(True)
#legend()
#show()
