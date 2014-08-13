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

errSH = [
0.321923,
0.219472,
0.105216,
0.0610161,
0.0269507,
0.0144433,
0.00730565,
0.00313838
]

errLX = [
0.13217,
0.0819829,
0.0354003,
0.0198031,
0.0171172,
0.0164807,
0.0162289,
0.0156743
]

errHeine1 = [
0.45,
0.39,
0.43,
0.48,
0.54,
0.59
]

errHeine2 = [
1,
0.59,
0.19,
0.063,
0.019,
0.0052,
0.0014
]

errHeine3 = [
1,
0.59,
0.19,
0.063,
0.019,
0.0052,
0.0014
]

errHeine4 = [
4.1,
0.36,
0.065,
0.0059,
0.0012,
0.00018,
2.80E-005
]


lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
plt.xlim([0.008,0.45]);

#plt.rc('font', family='serif')
#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSH, label=r"Weingarten", linewidth=4);
plt.loglog(hDEC, errLX, label=r"$\Delta\vec{x}$", linewidth=4);
plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM Degree 2", linewidth=4);
plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $H$ on ellipsoid");


grid(True)
legend(loc='lower right')
show()
