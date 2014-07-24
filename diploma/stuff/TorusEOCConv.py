from pylab import *
import matplotlib.pyplot as plt

eocL2 = [1.7472963007,1.8100340794,1.8517922781,1.9017333057,1.9356567635,1.957037723,1.9695867964,1.97929506];

eocMax = [1.676583294,1.7722152372,1.827141796,1.8875065925,1.9269600665,1.951624291,1.9656343129,1.9768600188]

h = [0.263248,0.214754,0.138422,0.0940568,0.0603108,0.0427842,0.0297978,0.0188374]


plt.plot(h, eocL2, "k-", label="L2", linewidth=3);
plt.plot(h, eocMax, "k--", label="Max", linewidth=3);

xlabel("h");
ylabel("EOC");
title("Torus");

grid(True)
legend()
show()
