from pylab import *
import matplotlib.pyplot as plt


hDEC = [
0.26332800,
0.21455300,
0.13811900,
0.09370030,
0.05996460,
0.04249460,
0.02957350,
0.01868440
]


errSK = [
4.2095500E-02,
2.7825800E-02,
1.1564200E-02,
5.4093000E-03,
2.2480700E-03,
1.1382300E-03,
5.5569800E-04,
2.2334900E-04
]

errGB = [
7.4500300E-03,
5.1209400E-03,
2.2474500E-03,
1.0716000E-03,
4.5101800E-04,
2.2980500E-04,
1.1287700E-04,
4.8314100E-05
]

errSKAvN = [
6.0133000E-02,
4.0723900E-02,
1.7549400E-02,
8.2933500E-03,
3.4722100E-03,
1.7641900E-03,
8.6257300E-04,
3.4712800E-04
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
#plt.xlim([0.007,0.8]);
plt.xlim([0.01,0.3]);

#plt.rc('font', family='serif')
#plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(hDEC, errSK, label=r"Weingarten", linewidth=4);
plt.loglog(hDEC, errSKAvN, label=r"Weingarten (AvN)", linewidth=4);
plt.loglog(hDEC, errGB, label=r"Gauss-Bonnet", linewidth=4);

xlabel("h");
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $K$ on torus");


grid(True)
legend(loc='lower right')
show()
