import csv
import sys
from pylab import *
import matplotlib.pyplot as plt

#fn = "output/sphereExtremeValues.csv"
fn = sys.argv[1]

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
        
        
scatter(t, phi)
grid(True)
#legend()
show()
