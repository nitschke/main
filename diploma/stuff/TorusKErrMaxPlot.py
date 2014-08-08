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
1.20032E-01,
8.15485E-02,
3.85582E-02,
1.76576E-02,
1.43730E-02,
1.00882E-02,
9.68054E-03,
1.07242E-02
]

errWeinConnN = [
1.13059E-01,
7.77378E-02,
3.62486E-02,
1.70546E-02,
7.49013E-03,
4.80814E-03,
4.22722E-03,
3.73022E-03
]

errWeinWN = [
1.18439E-01,
7.85310E-02,
3.63753E-02,
1.72349E-02,
7.43839E-03,
3.70765E-03,
1.88086E-03,
7.04270E-04
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
ylabel(r"$Err_{max}$");
title(r"$L_{\infty}$-Error for $K$ (and $\vec{\nu}$) on Torus");

grid(True)
legend(loc='lower right')
show()
