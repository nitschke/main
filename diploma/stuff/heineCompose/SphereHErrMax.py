from pylab import *
import matplotlib.pyplot as plt


hDEC = [
0.152342,
0.108972,
0.0663925,
0.0477329,
0.0305541,
0.0215183,
0.0152791,
0.00967015
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

errSH = [
2.852E-3,
1.471E-3,
5.492E-4,
2.843E-4,
1.166E-4,
5.787E-5,
2.917E-5,
1.169E-5
]

errLX = [
1.016E-5,
2.778E-6,
4.162E-7,
1.452E-7,
3.846E-8,
2.544E-8,
1.536E-8,
7.134E-9
]

errHeine1 = [
0.43,
0.47,
0.45,
0.44,
0.44,
0.48
]

errHeine2 = [
0.082,
0.031,
0.01,
0.0033,
0.001,
0.00033,
0.00011,
3.30E-005
]

errHeine3 = [
0.54,
0.23,
0.06,
0.02,
0.0061,
0.0019,
0.00061,
0.0002,
]

errHeine4 = [
0.26,
0.056,
0.0039,
0.0015,
0.00014,
3.80E-005,
7.50E-006,
2.70E-006
]


lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
plt.xlim([0.007,0.8]);


#plt.rc('font', family='serif')
#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSH, label=r"Weingarten", linewidth=4);
plt.loglog(hDEC, errLX, label=r"$\Delta\vec{x}$", linewidth=4);
plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM Degree 2", linewidth=4);
plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
#ylabel(r"$Err_{2}$");
#title(r"$L_{2}$-Error for $H$ on sphere");
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $H$ on sphere");

grid(True)
legend(loc='lower right')
show()
