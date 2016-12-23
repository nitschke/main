#!/usr/bin/python

import csv
from pylab import *
import sys



fnsDec = sys.argv[1:]

xlabel("time")
ylabel("zLoc")

for fn in fnsDec:
    timeDec = []
    zvalsDec = []
    
    with open(fn, 'rb') as f:
        reader = csv.reader(f, skipinitialspace = True)
        for row in reader:
            timeDec.append(row[0])
            zvalsDec.append(row[1])

    timeDec = array(timeDec, dtype=float)
    zvalsDec = array(zvalsDec, dtype=float)

    plot(timeDec, zvalsDec, "-o", label=fn, markevery=(timeDec.size-1))



legend(loc="upper right")

show()
