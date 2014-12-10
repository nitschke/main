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
4.934E-3,
2.527E-3,
9.389E-4,
4.855E-4,
1.990E-4,
9.869E-5,
4.975E-5,
1.993E-5
]

errGB = [
3.101E-3,
1.583E-3,
5.866E-4,
3.031E-4,
1.242E-4,
6.158E-5,
3.104E-5,
1.243E-5
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
0.0262689,
0.0165464,
0.00809494,
0.00498856,
0.00259324,
0.00155712,
0.000952146,
0.000500385
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

#plt.rc('xtick', labelsize=25)
#plt.rc('ytick', labelsize=25)
#plt.rc('axes', labelsize=25)
plt.xlim([0.007,0.16]);

plt.rcParams['font.size'] = 23
plt.rcParams['lines.linewidth'] = 6
plt.rcParams['lines.markersize'] = 12
plt.rcParams['lines.marker'] = 'd'

#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSK, label=r"W",);
plt.loglog(hDEC, errSKAvN, label=r"W,AvN");
plt.loglog(hDEC, errGB, label=r"GB");
#plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM-D2");
#plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
#plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $K$ on sphere");

grid(True, which="major", ls='-')
grid(True, which="minor")
legend(loc='lower right')
show()
