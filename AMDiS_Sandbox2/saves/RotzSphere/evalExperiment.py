#!/usr/bin/python

import csv
from pylab import *

hs = array([0.108972, 0.0663925, 0.0477329, 0.0305541, 0.0215183, 0.0152791])
ks = [2, 5, 10, 25, 50, 100]
fnsDec = []
fnsFem = []
for k in ks:
    fnsDec.append("dec/Tau01_" + str(k) + "k.RelKinErr.csv")
    fnsFem.append("fem/" + str(k) + "k/dissipationError.csv")

figure(1)
#subplot(211)
xlabel("time")
ylabel("Relative Kinetic Error")

finalRelKinErrDec = [];
for fn in fnsDec:
    timeDec = []
    RelKinErrDec = []
    
    with open(fn, 'rb') as f:
        reader = csv.reader(f, skipinitialspace = True)
        for row in reader:
            timeDec.append(row[0])
            RelKinErrDec.append(row[1])

    timeDec = array(timeDec, dtype=float)
    RelKinErrDec = array(RelKinErrDec, dtype=float)
    finalRelKinErrDec.append(RelKinErrDec[-1]);

    plot(timeDec, RelKinErrDec, label=fn)

gca().set_prop_cycle(None) #reset colors
finalRelKinErrFem = [];
for fn in fnsFem:
    timeFem = []
    RelKinErrFem = []

    with open(fn, 'rb') as f:
        reader = csv.reader(f, skipinitialspace = True)
        for row in reader:
            timeFem.append(row[0])
            RelKinErrFem.append(row[1])

    timeFem = array(timeFem, dtype=float)
    RelKinErrFem = array(RelKinErrFem, dtype=float)
    finalRelKinErrFem.append(RelKinErrFem[-1]);

    plot(timeFem, RelKinErrFem, "--", label=fn)


legend(loc="upper left")
#show()

#subplot(212)
figure(2)
xlabel("h_max")
ylabel("Relative Kinetic Error at t=10")

slopeDec, interceptDec = polyfit(log(hs), log(finalRelKinErrDec), 1)
loglog(hs, finalRelKinErrDec, "-*", label="final RelKinErr (DEC) (slope =" + str(slopeDec) + ")");

slopeFem, interceptFem = polyfit(log(hs), log(finalRelKinErrFem), 1)
loglog(hs, finalRelKinErrFem, "--*", label="final RelKinErr (FEM) (slope =" + str(slopeFem) + ")");


xs = hs
#loglog(xs, 1./sqrt(xs), label="slope 1/2")
#loglog(xs, xs, label="slope 1")
#loglog(xs, (xs)*(xs), label="slope 2")

legend(loc="lower right")
grid(True)
show()
