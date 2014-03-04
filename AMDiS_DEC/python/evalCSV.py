import csv
from pylab import *
import matplotlib.pyplot as plt

lineStyles = ['-','--', '-.', ':']
lw = 3;

#name = ["MaxAngleRatio", "AvAngleRatio"]
name = ["MaxDiameter", "AvDiameter"]
#name = ["MaxMaxAngle", "AvMaxAngle"]
#name = ["AvArea", "MinArea","MaxArea"]
n = len(name)
last = 10000
with open('../meshStatsSphereDivBy4.csv', 'rb') as f:
    reader = csv.DictReader(f)
    x = n*[ndarray((0,1),dtype=double)]
    k = 0
    for row in reader:
      for i in arange(n):
        x[i] = append(x[i], double(row[name[i]]))
      if k == last: break
      k = k+1


for i in arange(n):
  #semilogx(x[i], 'k'+lineStyles[i], label=name[i], linewidth=lw)
  plot(x[i], label=name[i])

#semilogx([1,last],[90,90],'k'+lineStyles[i+1],label="wellcentered", linewidth=lw)
#plot([0,last],[90,90],"--",label="wellcentered")

grid(True)
legend()
show()
