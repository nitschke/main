import csv
import sys
from pylab import *
import matplotlib.pyplot as plt

fn = "output.csv"

lw = 2;

name = ["t", "rhs", "phi", "phiExact"]
n = len(name)
with open(fn, 'rb') as f:
    reader = csv.DictReader(f)
    x = n*[ndarray((0,1),dtype=double)]
    for row in reader:
      for i in arange(n):
        x[i] = append(x[i], double(row[name[i]]))


for i in arange(1,n):
  plot(x[0], x[i], label=name[i], linewidth=lw)


grid(True)
legend()
show()

