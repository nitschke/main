from pylab import *
import matplotlib.pyplot as plt

err2dec = [
0.0152675,
0.00526172,
0.00363798,
0.00161156,
0.000772395,
0.000326629,
0.000166773,
0.000081775,
3.29872E-005
]

errmaxdec = [
0.0142554,
0.00513371,
0.0035787,
0.00160407,
0.00077352,
0.000328532,
0.0001681,
8.25597E-005,
3.33465E-005
]

err2fem = [
0.00421774,
0.000697127,
0.000424975,
0.000198093,
0.000109633,
5.25867E-005,
2.86589E-005,
1.47332E-005,
6.17951E-006
]

errmaxfem = [
0.00660004,
0.00112215,
0.000601834,
0.000303846,
0.000165767,
0.000076506,
4.07383E-005,
2.05796E-005,
0.000008506
]

h = [
0.484083,
0.263248,
0.214754,
0.138422,
0.0940568,
0.0603108,
0.0427842,
0.0297978,
0.0188374
]

plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.loglog(h, err2dec, label=r"$\mathrm{Err}_2$ (DEC)", linewidth=3);
plt.loglog(h, errmaxdec, label=r"$\mathrm{Err}_{\mathrm{max}}$ (DEC)", linewidth=3);
plt.loglog(h, err2fem, label=r"$\mathrm{Err}_2$ (FEM)", linewidth=3);
plt.loglog(h, errmaxfem, label=r"$\mathrm{Err}_{\mathrm{max}}$ (FEM)", linewidth=3);

xlabel("h");
ylabel("Err");
title("Torus");

grid(True)
legend(loc='best')
show()
