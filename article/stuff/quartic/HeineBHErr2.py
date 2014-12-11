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
0.0240724,
0.0120979
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
0.393089,
0.35968,
0.288263,
0.110784,
0.101324,
0.0845573,
0.0532833,
0.0265425,
0.00791621
]

errLX = [
5.2034900E-01,
4.6581700E-01,
3.8356700E-01,
8.9933000E-02,
4.6493600E-02,
3.2730100E-02,
1.8592300E-02,
8.9735700E-03,
2.6521100E-03
]

errHeine1 = [
0.17,
0.18,
0.18,
0.18,
0.18#,
#0.18
]

errHeine2 = [
0.021,
0.0055,
0.0014,
0.00034,
8.60E-005,
2.20E-005,
5.40E-006#,
#1.40E-006
]

errHeine3 = [
0.11,
0.048,
0.011,
0.0028,
0.0007,
0.00018,
4.50E-005#,
#1.10E-005
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
0.326261,
0.300841,
0.255107,
0.178652,
0.163547,
0.136016,
0.089295,
0.0472614,
0.014991
]


lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)


plt.xlim([0.01,0.5]);
plt.rcParams['font.size'] = 23
plt.rcParams['lines.linewidth'] = 6
plt.rcParams['lines.markersize'] = 12
plt.rcParams['lines.marker'] = 'd'

plt.loglog(hDEC, errSH, label=r"W");
plt.loglog(hDEC, errSHAvN, label=r"W,AvN");
plt.loglog(hDEC, errLX, label=r"LX");
#plt.loglog(hHeine1, errHeine1, label=r"FEM Degree 1", linewidth=4);
plt.loglog(hHeine2, errHeine2, label=r"FEM-D2", linewidth=4);
#plt.loglog(hHeine3, errHeine3, label=r"FEM Degree 3", linewidth=4);
#plt.loglog(hHeine4, errHeine4, label=r"FEM Degree 4", linewidth=4);

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $H$ on quartic surface");


grid(True, which="major", ls='-')
grid(True, which="minor")
legend(loc='lower right')
show()
