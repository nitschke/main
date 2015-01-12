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


hHeine2 = [
0.368,
0.204,
0.116,
0.0673,
0.0379,
0.0214
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


errHeine2 = [
0.25,
0.15,
0.073,
0.039,
0.018,
0.0085
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

hausdorff = [
3.2202400E-01,
2.4412700E-01,
1.3399600E-01,
2.4312500E-02,
1.4245800E-02,
7.8599600E-03,
3.2672400E-03,
1.2132800E-03
]

hausdorffTo = [
1.7253700E-01,
1.4136700E-01,
1.0111200E-01,
2.4732000E-02,
1.4233000E-02,
8.3140000E-03,
3.3690000E-03,
1.1360000E-03
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
plt.loglog(hHeine2, errHeine2, label=r"FEM-D2");
#plt.loglog(hDEC, hDEC, label=r"slope=1");
#plt.loglog(hDEC, array(hDEC)*array(hDEC), label=r"slope=2");
#plt.loglog(hDEC[0:-1], hausdorff, label=r"Hausdorff");
#plt.loglog(hDEC[0:-1], hausdorffTo, label=r"Hausdorff");
#plt.loglog(hausdorff, errSH[0:-1], label=r"W");
#plt.loglog(hausdorff, errSHAvN[0:-1], label=r"W,AvN");
#plt.loglog(hausdorff, errLX[0:-1], label=r"LX");

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $H$ on quartic surface");


grid(True, which="major", ls='-')
grid(True, which="minor")
legend(loc='lower right')
show()
