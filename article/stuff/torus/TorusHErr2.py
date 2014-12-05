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
1.5984600E-02,
1.1151400E-02,
4.9921500E-03,
2.3912400E-03,
1.0108400E-03,
5.1557800E-04,
2.5260900E-04,
1.0177600E-04
]

errLX = [
2.2616800E-03,
1.6409300E-03,
7.8803100E-04,
3.9930700E-04,
1.7499600E-04,
9.1056500E-05,
4.5227300E-05,
1.8688500E-05   
]

errSHAvN = [
2.7070000E-02,
1.8820000E-02,
8.4003400E-03,
4.0139300E-03,
1.6964100E-03,
8.6507200E-04,
4.2401000E-04,
1.7092600E-04
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
ylabel(r"$Err_{2}$");
title(r"$L_{2}$-Error for $H$ on torus");


grid(True)
legend(loc='lower right')
show()
