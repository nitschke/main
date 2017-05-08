from pylab import *
import scipy.sparse as sp
import scipy.sparse.linalg as la
import time

n=2 #modes
N=1000000 #grid points
p=10 # partititions

h = 1./(N-1)

fpsq= 4.*pi**2
hsq = h**2

f = lambda x: sin(2.*pi*n*x)
g = lambda x: fpsq*(n*n+1.)*f(x)

t1 = time.clock()
xh= zeros(p*N)
for j in arange(p):
    xh[j*N:(j+1)*N] = j + h*arange(N)
t2 = time.clock()

rhs = g(xh)
for j in arange(p):
    rhs[j*N] = rhs[j*N - 1] = 0.

#sm = zeros([p*N,p*N])
sm = sp.lil_matrix((p*N,p*N))
for j in arange(p):
    for i in arange(j*N+1,(j+1)*N-1):
        sm[i,i-1] = sm[i,i+1] = -1./hsq
        sm[i,i] = fpsq + 2./hsq

for j in arange(p):
    sm[j*N,j*N] = sm[j*N,j*N-1] = -3./2./h
    sm[j*N,j*N+1] = sm[j*N,j*N-2] = 4./2./h
    sm[j*N,j*N+2] = sm[j*N,j*N-3] = -1./2./h

for j in arange(p):
    sm[j*N-1,j*N] = 1.
    sm[j*N-1,j*N-1] = -1.
t3 = time.clock()
#smsp= sp.csr_matrix(sm)
t4 = time.clock()
print "NNZ: ", sm.nnz


#fhsol = solve(sm,rhs)
fhsol, info = la.bicgstab(sm,rhs,tol=1e-6, maxiter=10000)
t5 = time.clock()
print "Info: ", info

tauall = t5 - t1
print "mesh:  ", t2-t1, "s -> ", 100.*(t2-t1)/tauall, "%"
print "assem: ", t3-t2, "s -> ", 100.*(t3-t2)/tauall, "%"
print "conv:  ", t4-t3, "s -> ", 100.*(t4-t3)/tauall, "%"
print "solve: ", t5-t4, "s -> ", 100.*(t5-t4)/tauall, "%"

fh = f(xh)
plot(xh, fh, label="exact")
plot(xh, fhsol, label="solution")

legend()
#show()
