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


errSH = [
1.8696100E-02,
1.2918800E-02,
5.7317100E-03,
2.7275500E-03,
1.1522400E-03,
5.8760300E-04,
2.8819500E-04,
1.1608900E-04
]

errLX = [
4.5347900E-03,
3.5161400E-03,
1.7981200E-03,
9.3501400E-04,
4.1467400E-04,
2.1716400E-04,
1.0889600E-04,
5.6306700E-05
]

errSHAvN = [
3.8171600E-02,
2.6815100E-02,
1.1815700E-02,
5.7877400E-03,
2.3999000E-03,
1.2220400E-03,
5.9677600E-04,
2.4214400E-04
]

lineStyles = ['-','--', '-.', ':','', ' ', 'None']

plt.rc('text', usetex=True)

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
#plt.xlim([0.007,0.8]);
plt.xlim([0.01,0.3]);

#plt.rc('font', family='serif')
plt.loglog(hDEC, errSH, label=r"Weingarten", linewidth=4);
plt.loglog(hDEC, errSHAvN, label=r"Weingarten (AvN)", linewidth=4);
plt.loglog(hDEC, errLX, label=r"$\Delta\vec{x}$", linewidth=4);


xlabel("h");
ylabel(r"$Err_{\infty}$");
title(r"$L_{\infty}$-Error for $H$ on torus");


grid(True)
legend(loc='lower right')
show()
