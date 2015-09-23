import csv
from pylab import *
import matplotlib.pyplot as plt

lineStyles = ['-','--', '-.', ':']
lw = 3;

#name = ["MaxAngleRatio", "AvAngleRatio"]
#name = ["MaxDiameter", "AvDiameter"]
name = ["MaxMaxAngle", "AvMaxAngle"]
#name = ["AvArea", "MinArea","MaxArea"]
name2 = ["AvArea", "MinArea","MaxArea"]
n = len(name)
n2 = len(name2)
last = 100000
#with open('../meshStatsbunny.csv', 'rb') as f:
#with open('../meshStatsSphereDivBy4.csv', 'rb') as f:
fn = '../meshStatssphere0.6.csv'
with open(fn, 'rb') as f:
    reader = csv.DictReader(f)
    x = n*[ndarray((0,1),dtype=double)]
    x2 = n2*[ndarray((0,1),dtype=double)]
    k = 0
    for row in reader:
      for i in arange(n):
        x[i] = append(x[i], double(row[name[i]]))
      for i in arange(n2):
        x2[i] = append(x2[i], double(row[name2[i]]))
      if k == last: break
      k = k+1

fig, ax1 = plt.subplots()

for i in arange(n):
  #ax1.semilogx(x[i], label=name[i], linewidth=lw)
  ax1.semilogx(x[i], 'k'+lineStyles[i], label=name[i], linewidth=lw)
  #plot(x[i], label=name[i])

ax2 = ax1.twinx()
for i in arange(n2):
  ax2.semilogx(x2[i], label=name2[i], linewidth=lw)

#semilogx([1,last],[90,90],'k'+lineStyles[i+1],label="wellcentered", linewidth=lw)
#plot([0,last],[90,90],"--",label="wellcentered")

    
title(fn)
grid(True)
legend(loc='upper left')
show()
