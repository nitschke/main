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
1.23607,
1.51209,
0.938907,
0.480776,
0.432235,
0.434656,
0.304396
]

errGB = [
3.24383,
4.19967,
3.20761,
0.669981,
0.279146,
0.186159,
0.103139
]

errHeine1 = [
1.1,
1.2,
1.2,
1.2,
1.2,
1.2
]

errHeine2 = [
0.18,
0.13,
0.048,
0.016,
0.0053,
0.0017,
0.00053,
0.00017
]

errHeine3 = [
2.5,
0.52,
0.12,
0.037,
0.01,
0.0032,
0.001,
0.00034
]

errHeine4 = [
0.86,
0.14,
0.0081,
0.003,
0.00026,
7.50E-005,
1.50E-005,
5.40E-006
]

errSKAvN = [
0.838479,
1.01662,
0.681089,
0.645899,
0.637947,
0.572397,
0.423725
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
#ylabel(r"$Err_{2}$");
#title(r"$L_{2}$-Error for $K$ on quartic surface");
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $K$ on quartic surface");


grid(True)
legend(loc='lower right')
show()
