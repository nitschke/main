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

errSK = [
5.695E-3,
2.940E-3,
1.098E-3,
5.685E-4,
2.332E-4,
1.157E-4,
5.835E-5,
2.338E-5
]

errGB = [
3.583E-3,
1.844E-3,
6.871E-4,
3.555E-4,
1.458E-4,
7.233E-5,
3.647E-5,
1.461E-5
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


lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
plt.xlim([0.007,0.8]);


#plt.rc('font', family='serif')
#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSK, label=r"Weingarten", linewidth=4);
plt.loglog(hDEC, errGB, label=r"Gauss-Bonnet", linewidth=4);
plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM Degree 2", linewidth=4);
plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $K$ on sphere");

grid(True)
legend(loc='lower right')
show()
