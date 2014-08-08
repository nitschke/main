from pylab import *
import matplotlib.pyplot as plt

errNormal = [
1.12147E-02,
7.70817E-03,
3.37332E-03,
1.60779E-03,
7.11485E-04,
3.70873E-04,
1.89846E-04,
8.01637E-05
]


errWeinAvN = [
2.68038E-02,
1.86958E-02,
8.33410E-03,
4.06833E-03,
2.14116E-03,
2.01619E-03,
1.81599E-03,
1.29241E-03  
]

errWeinConnN = [
2.48099E-02,
1.73329E-02,
7.75109E-03,
3.74999E-03,
1.77822E-03,
1.40767E-03,
1.18653E-03,
8.23916E-04
]

errWeinWN = [
1.59035E-02,
1.11176E-02,
4.96833E-03,
2.38374E-03,
1.00729E-03,
5.13812E-04,
2.51673E-04,
1.01450E-04
]

h = [
0.263248,
0.214754,
0.138422,
0.0940568,
0.0603108,
0.0427842,
0.0297978,
0.0188374
]

lineStyles = ['-','--', '-.', ':']

plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.loglog(h, errWeinWN, 'k-', label=r"(HExN)", linewidth=4);
plt.loglog(h, errWeinAvN, 'k--', label=r"(HAvN)", linewidth=4);
plt.loglog(h, errNormal, 'k:', label=r"(AvN)", linewidth=4);
plt.loglog(h, errWeinConnN, 'k-.', label=r"(HConnN)", linewidth=4);

xlabel("h");
ylabel(r"$Err_2$");
title(r"$L_2$-Error for $H$ (and $\vec{\nu}$) on Torus");

grid(True)
legend(loc='lower right')
show()
