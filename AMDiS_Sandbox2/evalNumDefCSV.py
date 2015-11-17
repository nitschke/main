import csv
import sys
from pylab import *
import matplotlib.pyplot as plt

#fn = "output/sphereEnergies.csv"
fn = sys.argv[1]

lineStyles = ['-','--', '-.', ':']
lw = 3;

name = ["Time", "NumberOfDefects"]
n = len(name)
last = 1000000000
breakAt = 2
with open(fn, 'rb') as f:
    reader = csv.DictReader(f, skipinitialspace = True)
    x = n*[ndarray((0,1),dtype=double)]
    k = 0
    for row in reader:
      for i in arange(n):
        x[i] = append(x[i], double(row[name[i]]))
      #if int(row["NumberOfDefects"]) <= breakAt: break
      if k == last: break
      k = k+1


for i in arange(1,n):
  #semilogx(x[0], x[i], 'k'+lineStyles[i], label=name[i], linewidth=lw)
  #plot(x[0], x[i], 'k'+lineStyles[i], label=name[i], linewidth=lw)
  #semilogy(x[0], x[i], label=name[i], linewidth=lw)
  #semilogx(x[0], x[i], label=name[i], linewidth=lw)
  loglog(x[0], x[i], label=name[i], linewidth=lw)
  #plot(x[0], x[i], label=name[i], linewidth=lw)
    

grid(True)
legend()
show()
