import csv
import sys
from pylab import *
import matplotlib.pyplot as plt

with open('virus_kS_intK2_vol_length.csv', 'rb') as f:
    reader = csv.DictReader(f, skipinitialspace = True)
    kS = ndarray((0,1),dtype=double)
    intK2 = ndarray((0,1),dtype=double)
    vol = ndarray((0,1),dtype=double)
    length = ndarray((0,1),dtype=double)
    for row in reader:
        kS = append(kS, double(row['kS']))
        intK2 = append(intK2, double(row['intK2']))
        vol = append(vol, double(row['vol']))
        length = append(length, double(row['length']))

f, axarr = plt.subplots(4, sharex=True)

axarr[0].plot(kS,intK2)
axarr[0].set_title('int(K^2)')
axarr[1].plot(kS,vol)
axarr[1].set_title('vol')
axarr[2].plot(kS,length)
axarr[2].set_title('length')
axarr[3].plot(kS,intK2/vol)
axarr[3].set_title('intK2/vol')

plt.xlabel('kS')

plt.show()
