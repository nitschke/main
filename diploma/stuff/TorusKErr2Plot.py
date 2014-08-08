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
6.00992E-02,
4.06086E-02,
1.77129E-02,
8.65477E-03,
4.12595E-03,
2.78750E-03,
2.18022E-03,
1.60895E-03
]

errWeinConnN = [
5.89887E-02,
4.01822E-02,
1.75346E-02,
8.38940E-03,
3.64332E-03,
2.12024E-03,
1.39385E-03,
8.84444E-04
]

errWeinWN = [
4.21395E-02,
2.76862E-02,
1.16123E-02,
5.42485E-03,
2.25709E-03,
1.14510E-03,
5.61061E-04,
2.26305E-04
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
plt.loglog(h, errWeinWN, 'k-', label=r"(KExN)", linewidth=4);
plt.loglog(h, errWeinAvN, 'k--', label=r"(KAvN)", linewidth=4);
plt.loglog(h, errNormal, 'k:', label=r"(AvN)", linewidth=4);
plt.loglog(h, errWeinConnN, 'k-.', label=r"(KConnN)", linewidth=4);

xlabel("h");
ylabel(r"$Err_2$");
title(r"$L_2$-Error for $K$ (and $\vec{\nu}$) on Torus");

grid(True)
legend(loc='lower right')
show()
