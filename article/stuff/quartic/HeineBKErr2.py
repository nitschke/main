from pylab import *
import matplotlib.pyplot as plt


hDEC = [
0.413931,
0.339603,
0.257294,
0.109466,
0.078606,
0.0598497,
0.0390528,
0.02407240,
0.01209790
]

hHeine2 = [
0.368,
0.204,
0.116,
0.0673,
0.0379,
0.0214
]


errSK = [
9.8203400E-01,
1.0218700E+00,
8.8338600E-01,
2.2116600E-01,
2.2545500E-01,
2.0954000E-01,
1.3965000E-01,
7.1082000E-02,
2.2288900E-02
]

errGB = [
1.9160700E+00,
1.8700100E+00,
1.8865600E+00,
4.4229600E-01,
1.7534500E-01,
1.0438100E-01,
5.1754700E-02,
2.4303200E-02,
7.2582200E-03
]


errHeine2 = [
2.7,
4.4,
2.4,
0.96,
0.34,
0.11
]


errSKAvN = [
7.6599900E-01,
7.5564000E-01,
6.2178100E-01,
3.7265700E-01,
3.7807900E-01,
3.2502400E-01,
2.2200500E-01,
1.2201300E-01,
4.0263700E-02
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.xlim([0.01,0.5]);

plt.rcParams['font.size'] = 23
plt.rcParams['lines.linewidth'] = 6
plt.rcParams['lines.markersize'] = 12
plt.rcParams['lines.marker'] = 'd'

plt.loglog(hDEC, errSK, label=r"W");
plt.loglog(hDEC, errSKAvN, label=r"W,AvN");
plt.loglog(hDEC, errGB, label=r"GB");
plt.loglog(hHeine2, errHeine2, label=r"FEM-D2");

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $K$ on quartic surface");


grid(True, which="major", ls='-')
grid(True, which="minor")
legend(loc='lower right')
show()
