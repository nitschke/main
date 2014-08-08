from pylab import *
import matplotlib.pyplot as plt

errNormal = [
1.34256E-02,
9.38635E-03,
4.04736E-03,
2.10221E-03,
1.22951E-03,
6.88432E-04,
4.27560E-04,
2.22600E-04
]


errWeinAvN = [
3.77115E-02,
2.74883E-02,
1.22584E-02,
7.67575E-03,
8.90489E-03,
6.35418E-03,
6.07699E-03,
5.49763E-03
]

errWeinConnN = [
2.92189E-02,
2.11396E-02,
9.44613E-03,
5.72605E-03,
5.85822E-03,
3.94152E-03,
3.31123E-03,
2.51031E-03
]

errWeinWN = [
1.83957E-02,
1.27709E-02,
5.64408E-03,
2.71495E-03,
1.15397E-03,
5.90779E-04,
2.93789E-04,
1.24545E-04
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
ylabel(r"$Err_{max}$");
title(r"$L_{\infty}$-Error for $H$ (and $\vec{\nu}$) on Torus");

grid(True)
legend(loc='lower right')
show()
