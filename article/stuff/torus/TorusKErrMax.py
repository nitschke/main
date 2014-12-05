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
1.1558300E-01,
8.2688500E-02,
3.5303500E-02,
1.7127800E-02,
7.1558800E-03,
3.7181200E-03,
1.8430000E-03,
7.0564400E-04
]

errGB = [
1.1173400E-02,
7.9415900E-03,
3.6171700E-03,
1.7750600E-03,
7.5732700E-04,
3.8986400E-04,
2.0112400E-04,
1.2669400E-04
]

errSKAvN = [
1.1867800E-01,
8.2649900E-02,
3.8006400E-02,
1.7982600E-02,
7.7854000E-03,
3.7474100E-03,
1.9245100E-03,
7.1923200E-04
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
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $K$ on torus");


grid(True)
legend(loc='lower right')
show()
