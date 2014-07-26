from pylab import *
import matplotlib.pyplot as plt

err2dec = [
0.00451724,
0.00230461,
0.00085386,
0.00044110,
0.00018068,
0.00008960,
0.00004517,
0.00001809
]

errmaxdec = [
0.00467138,
0.00237221,
0.00088076,
0.00045522,
0.00018645,
0.00009249,
0.00004662,
0.00001868
]

err2fem = [
0.00701125,
0.00358810,
0.00133213,
0.00068858,
0.00028215,
0.00013995,
0.00007055,
0.00002826
]

errmaxfem = [
0.00865265,
0.00455121,
0.00176727,
0.00093964,
0.00039924,
0.00020360,
0.00010537,
0.00004367
]

h = [
0.152342,
0.108972,
0.0663925,
0.0477329,
0.0305541,
0.0215183,
0.0152791,
0.00967015
]

plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.loglog(h, err2dec, label=r"$\mathrm{Err}_2$ (DEC)", linewidth=3);
plt.loglog(h, errmaxdec, label=r"$\mathrm{Err}_{\mathrm{max}}$ (DEC)", linewidth=3);
plt.loglog(h, err2fem, label=r"$\mathrm{Err}_2$ (FEM)", linewidth=3);
plt.loglog(h, errmaxfem, label=r"$\mathrm{Err}_{\mathrm{max}}$ (FEM)", linewidth=3);

xlabel("h");
ylabel("Err");
title("Sphere");

grid(True)
legend(loc='best')
show()
