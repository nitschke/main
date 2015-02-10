from pylab import *

n = 100

r = rand(n)

A = zeros([n, 2])

A[:,0] = r*r
A[:,1] = r

b = ones_like(r)

c = inv(A.T.dot(A)).dot(A.T.dot(b))

def poly(x):
    return c[0]*x*x + c[1]*x - 1.0

plot(r, zeros_like(r), "*")

x = linspace(r.min(), r.max(), 100)
plot(x, poly(x))

show()
