from pylab import *
import matplotlib.pyplot as plt

eocL2 = [2.0092942158,
2.0040930938,
2.0018345526,
2.0008089195,
2.0003849111,
2.000437711,
1.9999448405];

eocMax = [2.0225601077,
1.9995573186,
2.0002292995,
2.0008136244,
1.999641636,
2.0004803364,
1.9997695362,
];

h = [0.108972,
0.0663925,
0.0477329,
0.0305541,
0.0215183,
0.0152791,
0.00967015];


plt.plot(h, eocL2, "k-", label="L2", linewidth=3);
plt.plot(h, eocMax, "k--", label="Max", linewidth=3);

xlabel("h");
ylabel("EOC");
title("Sphere");

grid(True)
legend()
show()
