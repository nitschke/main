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
0.0968374,
0.0639814,
0.0307902,
0.0175462,
0.00775927,
0.00415056,
0.00209925,
0.00089715
]

errLX = [
0.037736,
0.0224912,
0.0103995,
0.00589007,
0.0026113,
0.00140573,
0.000726692,
0.00032242
]

errHeine1 = [
0.19,
0.18,
0.18,
0.18,
0.18,
0.18
]

errHeine2 = [
0.23,
0.09,
0.025,
0.0067,
0.0017,
0.00044,
0.00011
]

errHeine3 = [
0.23,
0.09,
0.025,
0.0067,
0.0017,
0.00044,
0.00011
]

errHeine4 = [
0.39,
0.04,
0.0056,
0.00063,
8.00E-005,
1.00E-005,
1.30E-006
]

errSHAvN = [
0.160546,
0.111961,
0.0576774,
0.0341154,
0.0156266,
0.00849949,
0.00435455,
0.00188527
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
plt.loglog(hDEC, errSHAvN, label=r"Weingarten (AvN)", linewidth=4);
plt.loglog(hDEC, errLX, label=r"$\Delta\vec{x}$", linewidth=4);
plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM Degree 2", linewidth=4);
plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $H$ on ellipsoid");


grid(True)
legend(loc='lower right')
show()
