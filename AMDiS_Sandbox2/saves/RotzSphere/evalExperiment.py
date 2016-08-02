#!/usr/bin/python

import csv
from pylab import *

ks = [2, 5, 10, 25, 50]
fns = []
for k in ks:
    fns.append("Tau01_" + str(k) + "k.RelKinErr.csv")

figure(1)
#subplot(211)
xlabel("time")
ylabel("Relative Kinetic Error")

finalRelKinErr = [];
for fn in fns:
    time = []
    RelKinErr = []
    
    with open(fn, 'rb') as f:
        reader = csv.reader(f, skipinitialspace = True)
        for row in reader:
            time.append(row[0])
            RelKinErr.append(row[1])

    time = array(time, dtype=float)
    RelKinErr = array(RelKinErr, dtype=float)
    finalRelKinErr.append(RelKinErr[-1]);

    plot(time, RelKinErr, label=fn)


legend(loc="upper left")
#show()

#subplot(212)
figure(2)
xlabel("Meshidentifier")
ylabel("Relative Kinetic Error at t=10")

loglog(ks, finalRelKinErr, "-*", label="final RelKinErr");
xs = arange(2,100,1)
loglog(xs, 1./sqrt(xs), label="slope 1/2")
loglog(xs, 1./xs, label="slope 1")
loglog(xs, 1./(xs*xs), label="slope 2")

legend()
grid(True)
show()
