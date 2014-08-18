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
2.470E-3,
1.264E-3,
4.696E-4,
2.428E-4,
9.948E-5,
4.935E-5,
2.488E-5,
9.966E-6
]

errLX = [
4.550E-6,
1.352E-6,
2.108E-7,
6.034E-8,
1.217E-8,
3.887E-9,
1.549E-9,
5.137E-10
]

errHeine1 = [
0.17,
0.18,
0.18,
0.18,
0.18,
0.18
]

errHeine2 = [
0.021,
0.0055,
0.0014,
0.00034,
8.60E-005,
2.20E-005,
5.40E-006,
1.40E-006
]

errHeine3 = [
0.11,
0.048,
0.011,
0.0028,
0.0007,
0.00018,
4.50E-005,
1.10E-005
]

errHeine4 = [
0.058,
0.012,
0.0008,
0.00014,
9.70E-006,
1.40E-006,
1.20E-007,
1.70E-008
]

errSHAvN = [
0.0132212,
0.00830892,
0.00405554,
0.00249747,
0.00129781,
0.000779469,
0.000476943,
0.000250996
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
plt.loglog(hDEC, errSHAvN, label=r"Weingarten (AvN)", linewidth=4);
plt.loglog(hDEC, errLX, label=r"$\Delta\vec{x}$", linewidth=4);
plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM Degree 2", linewidth=4);
plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $H$ on sphere");

grid(True)
legend(loc='lower right')
show()
