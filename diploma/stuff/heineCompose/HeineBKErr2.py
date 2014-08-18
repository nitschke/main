from pylab import *
import matplotlib.pyplot as plt


hDEC = [
0.413931,
0.339603,
0.257294,
0.109466,
0.078606,
0.0598497,
0.0390528
]

hHeine1 = [
0.126,
0.0717,
0.0407,
0.0231,
0.0131,
0.00743
]

hHeine2 = [
0.395,
0.242,
0.138,
0.0784,
0.0441,
0.0249,
0.014,
0.00787
]

hHeine3 = [
0.515,
0.284,
0.154,
0.0863,
0.0486,
0.0273,
0.0154,
0.00856
]

hHeine4 = [
0.793,
0.448,
0.244,
0.144,
0.0817,
0.0464,
0.0261,
0.0146
]

errSK = [
0.982034,
1.02187,
0.883386,
0.221166,
0.225455,
0.20954,
0.13965
]

errGB = [
1.91607,
1.87001,
1.88656,
0.442296,
0.175345,
0.104381,
0.0517547
]

errHeine1 = [
0.41,
0.41,
0.42,
0.42,
0.42,
0.42
]

errHeine2 = [
0.052,
0.025,
0.0073,
0.0019,
0.00049,
0.00012,
3.10E-005,
7.70E-006
]

errHeine3 = [
0.27,
0.096,
0.021,
0.0048,
0.0012,
0.00031,
7.80E-005,
2.00E-005
]

errHeine4 = [
0.14,
0.025,
0.0015,
0.00028,
1.90E-005,
2.80E-006,
2.40E-007,
3.40E-008
]

errSKAvN = [
0.765999,
0.75564,
0.621781,
0.372657,
0.378079,
0.325024,
0.222005
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
plt.xlim([0.007,0.8]);

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
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $K$ on quartic surface");


grid(True)
legend(loc='lower right')
show()
