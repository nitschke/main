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
0.0159,
0.00833
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
0.47307,
0.326632,
0.157128,
0.0915163,
0.0404486,
0.0216735,
0.0109678,
0.00471233
]

errGB = [
0.0585884,
0.0552849,
0.053577,
0.0517027,
0.0512228,
0.052223,
0.0519484,
0.0534929
]

errHeine1 = [
1.4,
1.7,
1.7,
1.7,
1.7,
1.8
]

errHeine2 = [
8.3,
3.4,
1.2,
0.31,
0.086,
0.024,
0.0063,
]

errHeine3 = [
8.3,
3.4,
1.2,
0.31,
0.086,
0.024,
0.0063
]

errHeine4 = [
72,
2.5,
0.24,
0.027,
0.004,
0.00057,
9.90E-005
]

errSKAvN = [
0.603145,
0.464826,
0.262582,
0.166081,
0.0790886,
0.0625738,
0.0611925,
0.0593363
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
plt.xlim([0.008,0.45]);

#plt.rc('font', family='serif')
#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSK, label=r"Weingarten", linewidth=4);
plt.loglog(hDEC, errSKAvN, label=r"Weingarten (AvN)", linewidth=4);
plt.loglog(hDEC, errGB, label=r"Gauss-Bonnet", linewidth=4);
plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM Degree 2", linewidth=4);
plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $K$ on ellipsoid");


grid(True)
legend(loc='lower right')
show()
