from pylab import *
import matplotlib.pyplot as plt

errNormal2 = [

]

errNormalMax = [

]

errWein2 = [

]

errWeinMax = [

]

errWeinWN2 = [

]

errWeinWNMax = [

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
plt.loglog(h, errWein, 'k-', label=r"$\mathrm{Err}_2$ (DEC)", linewidth=4);
plt.loglog(h, errWeinWN, 'k--', label=r"$\mathrm{Err}_{\mathrm{max}}$ (DEC)", linewidth=4);
plt.loglog(h, err, 'k-.', label=r"$\mathrm{Err}_2$ (FEM)", linewidth=4);
plt.loglog(h, errmaxfem, 'k:', label=r"$\mathrm{Err}_{\mathrm{max}}$ (FEM)", linewidth=4);

xlabel("h");
ylabel("Err");
title("Torus");

grid(True)
legend(loc='best')
show()
