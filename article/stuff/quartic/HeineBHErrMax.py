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

hHeine1 = [
0.126,
0.0717,
0.0407,
0.0231,
0.0131#,
#0.00743
]

hHeine2 = [
0.395,
0.242,
0.138,
0.0784,
0.0441,
0.0249,
0.014#,
#0.00787
]

hHeine3 = [
0.515,
0.284,
0.154,
0.0863,
0.0486,
0.0273,
0.0154#,
#0.00856
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
5.0919400E-01,
4.5445000E-01,
3.7389200E-01,
3.5854100E-01,
3.7374000E-01,
3.2585500E-01,
2.2238800E-01,
1.1347300E-01,
3.3842200E-02
]

errLX = [
7.9905700E-01,
6.6901300E-01,
6.5007700E-01,
2.3083600E-01,
1.5406000E-01,
1.3104200E-01,
8.3735300E-02,
3.9446400E-02,
1.1085200E-02
]

errHeine1 = [
0.43,
0.47,
0.45,
0.44,
0.44#,
#0.48
]

errHeine2 = [
0.082,
0.031,
0.01,
0.0033,
0.001,
0.00033,
0.00011#,
#3.30E-005
]

errHeine3 = [
0.54,
0.23,
0.06,
0.02,
0.0061,
0.0019,
0.00061#,
#0.0002
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

errSHAvN = [
4.1263900E-01,
4.4516200E-01,
5.0383600E-01,
5.0040200E-01,
5.1129600E-01,
4.4588000E-01,
3.2250500E-01,
1.8432000E-01,
6.2137800E-02
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
#plt.xlim([0.007,0.8]);
plt.xlim([0.01,0.8]);

#plt.rc('font', family='serif')
#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSH, label=r"Weingarten", linewidth=4);
plt.loglog(hDEC, errSHAvN, label=r"Weingarten (AvN)", linewidth=4);
#plt.loglog(hDEC, errGB, label=r"Gauss-Bonnet", linewidth=4);
plt.loglog(hDEC, errLX, label=r"$\Delta\vec{x}$", linewidth=4);
plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM Degree 2", linewidth=4);
plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $H$ on quartic surface");


grid(True)
legend(loc='lower right')
show()
