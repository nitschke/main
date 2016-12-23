#!/usr/bin/python

import csv
from pylab import *
import sys


fnsDec = sys.argv[1:]

try:
    maxtime = float(fnsDec[0])
except ValueError:
    maxtime = inf

if maxtime < inf:
    print maxtime
    fnsDec = fnsDec[1:]



xlabel("time")
ylabel("kinetic error")

for fn in fnsDec:
    timeDec = []
    KinErrDec = []
    
    with open(fn, 'rb') as f:
        reader = csv.reader(f, skipinitialspace = True)
        for row in reader:
            if float(row[0]) > maxtime:
                break
            timeDec.append(row[0])
            KinErrDec.append(row[1])

    timeDec = array(timeDec, dtype=float)
    KinErrDec = array(KinErrDec, dtype=float)

    plot(timeDec, KinErrDec, "-o", label=fn, markevery=(timeDec.size-1))



legend(loc="upper right")

show()
