from triangle import *
from pylab import *

a = array([ 0 ,  0 , -1 ])
b = array([ 0.707107 , 0 , -0.707107 ])
c = array([ 0.325058 , -0.325058 , -0.888074])

t = triangle(a, b ,c)
print t
t.plot()
