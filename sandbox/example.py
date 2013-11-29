from triangle import *
from pylab import *

a = array([ 0.,  0.,  0.])
b = array([ 1.,  2.,  3.])
c = array([ 1.,  3.,  1.])

t = triangle(a, b ,c)
print t
t.plot()
