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

errSK = [
9.8203400E-01,
1.0218700E+00,
8.8338600E-01,
2.2116600E-01,
2.2545500E-01,
2.0954000E-01,
1.3965000E-01,
7.1082000E-02,
2.2288900E-02
]

errGB = [
1.9160700E+00,
1.8700100E+00,
1.8865600E+00,
4.4229600E-01,
1.7534500E-01,
1.0438100E-01,
5.1754700E-02,
2.4303200E-02,
7.2582200E-03
]

errHeine1 = [
0.41,
0.41,
0.42,
0.42,
0.42#,
#0.42
]

errHeine2 = [
0.052,
0.025,
0.0073,
0.0019,
0.00049,
0.00012,
3.10E-005#,
#7.70E-006
]

errHeine3 = [
0.27,
0.096,
0.021,
0.0048,
0.0012,
0.00031,
7.80E-005#,
#2.00E-005
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
7.6599900E-01,
7.5564000E-01,
6.2178100E-01,
3.7265700E-01,
3.7807900E-01,
3.2502400E-01,
2.2200500E-01,
1.2201300E-01,
4.0263700E-02
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
