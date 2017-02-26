from pylab import *
import mpmath

mpmath.mp.dps = 50     # higher precision for demonstration
#Nv=100000
Nv=100000
a = [mpmath.sin(2.*mpmath.pi*n/Nv) for n in range(Nv)]
b = array(a)
print  b.dot(b) - Nv/2

b32= array(a,dtype=float32)
b64= array(a,dtype=float64)
b128= array(a,dtype=float128)

print  b32.dot(b32) - Nv/2.
print  b64.dot(b64) - Nv/2.
print  b128.dot(b128) - Nv/2.

print
print mpmath.norm(b) - mpmath.sqrt(Nv/2)
print norm(b32)  - sqrt(Nv/2)
print norm(b64)  - sqrt(Nv/2)
print norm(b128)  -sqrt(Nv/2)
