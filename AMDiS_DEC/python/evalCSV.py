import csv
from pylab import *
#name = ["MaxMaxAngle", "AvMaxAngle"]
name = ["AvArea", "MinArea","MaxArea"]
n = len(name)
last = 20
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
  #semilogx(x[i], label=name[i])
  plot(x[i], label=name[i])

legend()
show()
