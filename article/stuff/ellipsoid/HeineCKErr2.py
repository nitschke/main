from pylab import *
import matplotlib.pyplot as plt


hDEC = [
0.249117,
0.178909,
0.109365,
0.0784945,
0.0503854,
0.0364551,
0.0257207,
0.0167543
]

hHeine1 = [
0.2,
0.117,
0.0674,
0.0384,
0.0218,
0.0127
]

hHeine2 = [
0.295,
0.19,
0.107,
0.0588,
0.0309,
0.0159
]

hHeine3 = [
0.309,
0.198,
0.131,
0.0788,
0.0457,
0.0259,
0.0145
]

hHeine4 = [
0.443,
0.265,
0.153,
0.0908,
0.0527,
0.0297,
0.0168
]

errSK = [
0.209717,
0.139978,
0.0689603,
0.0397115,
0.0177063,
0.00952316,
0.00482189,
0.00206855
]

errGB = [
0.0605908,
0.0304249,
0.0117452,
0.00648864,
0.00307695,
0.00183794,
0.00111634,
0.0006237
]

errHeine1 = [
0.76,
0.83,
0.83,
0.83,
0.83,
0.84
]

errHeine2 = [
0.69,
0.49,
0.21,
0.069,
0.025,
0.011
]

errHeine3 = [
2.2,
0.76,
0.16,
0.041,
0.011,
0.0027,
0.00067
]

errHeine4 = [
7.5,
0.23,
0.024,
0.0027,
0.00034,
4.30E-005,
5.50E-006
]

errSKAvN = [
0.321983,
0.232983,
0.125849,
0.0759988,
0.0353994,
0.0194333,
0.00999532,
0.00436155
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.xlim([0.015,0.3]);
plt.rcParams['font.size'] = 23
plt.rcParams['lines.linewidth'] = 6
plt.rcParams['lines.markersize'] = 12
plt.rcParams['lines.marker'] = 'd'

#plt.rc('font', family='serif')
#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSK, label=r"W");
plt.loglog(hDEC, errSKAvN, label=r"W,AvN");
plt.loglog(hDEC, errGB, label=r"GB");
#plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM-D2");
#plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
#plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $K$ on ellipsoid");


grid(True, which="major", ls='-')
grid(True, which="minor")
legend(loc='lower right')
show()
