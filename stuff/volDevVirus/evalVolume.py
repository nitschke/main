from pylab import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

kS = [10,15,20,25,30,35,40,45,50,55,60,65,70,100,200]
volPi = [4.02426,3.9347,3.8615,3.79879,3.74362,3.6942,3.64932,3.60813,3.57001,3.53452,3.50129,3.47002,3.44047,3.29002,2.94474]

f = UnivariateSpline(kS, volPi, k=4)
kSNew=np.linspace(kS[1], kS[-1], num=100, endpoint=True)

plt.plot(kS,volPi,'o')
plt.plot(kS,volPi)

#plt.plot(kSNew,f(kSNew))


plt.show()
