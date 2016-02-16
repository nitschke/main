import csv
from pylab import *
import matplotlib.pyplot as plt

## Old Data (coarser time discretation)
#Data = [
#        [0.5, 22.95678499, 0.21],
#        [0.6, 22.89465497, 0.35],
#        [0.7, 22.86525998, 0.55],
#        [0.8, 22.87307081, 0.88],
#        [0.9, 22.92027478, 1.52],
#        [1.0, 23.00728656, 3.43],
#        [1.01, 23.01790565, 4.39],
#        [1.02, 23.02892059, 5.3],
#        [1.03, 23.04032154, 5.38],
#        [1.04, 23.05197395, 6.32],
#        [1.05, 23.06419045, 7.30],
#        [1.06, 23.0766847, 8.32],
#        [1.07, 23.0894524, 10.31],
#        [1.08, 23.10262915, 13.44],
#        [1.09, 23.11619359, 24.43],
#        [1.1, 33.00051755, -1],
#        [1.2, 32.85992495, -1],
#        [1.3, 32.70922517, -1],
#        [1.4, 32.57215187, -1],
#        [1.5, 32.40689058, -1]
#       ]

# [C, final energy, fusion time]
Data = [
        [0.5, 22.95678393, 0.21],
        [0.6, 22.8946556, 0.35],
        [0.7, 22.86526009, 0.56],
        [0.8, 22.87307092, 0.89],
        [0.9, 22.92027424, 1.55],
        [1.0, 23.00728655, 3.29],
        [1.01, 23.01790565, 3.64],
        [1.02, 23.02891986, 4.07],
        [1.03, 23.04032154, 4.60],
        [1.04, 23.05197434, 5.30],
        [1.05, 23.06419045, 6.18],
        [1.06, 23.07668505, 7.80],
        [1.07, 23.08945241, 10.24],
        [1.08, 23.10262949, 13.22],
        [1.09, 23.11619357, 24.48],
        [1.091, 23.11770407, 27.56],
        [1.092, 23.11904568, 31.44],
        [1.093, 23.12038767, 37.79],
        [1.094, 23.12173266, 48.89],
        [1.095, 23.12308115, 85.49],
        [1.096, 33.00629998, -1],
        [1.097, 33.00556329, -1],
        [1.098, 33.00430197, -1],
        [1.099, 33.00296753, -1],
        [1.1, 33.00051786, -1],
        [1.2, 32.85992495, -1],
        [1.3, 32.70922511, -1],
        [1.4, 32.57215184, -1],
        [1.5, 32.40689058, -1]
    ]


DataT = array(Data).transpose()
plt.plot(DataT[0],DataT[1], linestyle="-", marker="d", linewidth=3, label="Final Energy E")

DataFilteredT = array(filter(lambda dat: dat[2] > 0, Data)).transpose()
plt.plot(DataFilteredT[0],DataFilteredT[2], linestyle="-", marker="o", linewidth=3, label="Fusion Time t")

plt.axvspan(0.5, 1.095, facecolor='0.5', alpha=0.5)
plt.text(0.8, 10, "2 Defects", horizontalalignment='center')
plt.text(1.3, 10, "4 Defects", horizontalalignment='center')

xlim(0.5,1.5)
xlabel("Shift C")
ylabel("Energy E / Time t")
title("Nonic Surfaces (with fixed proportion factor r=0.95)")
legend()
plt.show()