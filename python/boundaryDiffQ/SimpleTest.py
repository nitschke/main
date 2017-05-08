from pylab import *

n=2 #modes
N=20 #grid points

h = 1./(N-1)

fpsq= 4.*pi**2
hsq = h**2

f = lambda x: sin(2.*pi*n*x)
g = lambda x: fpsq*(n*n+1.)*f(x)

fh = f(h*arange(N))

rhs = g(h*arange(N))
rhs[-1] = 0.0 #coupling condition

sm = zeros([N,N])
sm[0,0] = fpsq - 2./hsq
sm[0,1] = 5./hsq
sm[0,2] = -4./hsq
sm[0,3] = 1./hsq
for i in arange(1,N-1):
    sm[i,i-1] = sm[i,i+1] = -1./hsq
    sm[i,i] = fpsq + 2./hsq
sm[-1,0] = 1.
sm[-1,-1] = -1

fhsol = solve(sm,rhs)

smcl= sm
smcl[0,0] = fpsq + 2./hsq
smcl[0,-1] = smcl[0,1] = -1./hsq
smcl[0,2] = smcl[0,3] = 0.
#fhsolcl = solve(smcl[:-1,:-1],rhs[:-1])
fhsolcl = solve(smcl,rhs)

rhsed= rhs
rhsed[0] = 0.
smed = sm
smed[0,0] = smed[0,-1] = -3./2./h
smed[0,1] = smed[0,-2] = 4./2./h
smed[0,2] = smed[0,-3] = -1./2./h
smed[0,3] = 0.
fhsoled= solve(smed,rhsed)

xh = h*arange(N)
plot(xh, fh, label="exact")
plot(xh, fhsol, label="solution")
plot(xh, fhsoled, label="solution eq deriv")
#plot(xh[:-1], fhsolcl, label="solution classic")
plot(xh, fhsolcl, label="solution classic")

legend()
show()
