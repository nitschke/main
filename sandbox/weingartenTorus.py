from pylab import *

def II(t,p):
    R = 2.0
    r = 0.5
    x_t = (R+r*cos(p)) * matrix([-sin(t), cos(t), 0.]).T
    x_p = (-r) * matrix([cos(t)*sin(p), sin(t)*sin(p), -cos(p)]).T
    n_t = (1./((R+r*cos(p))*(R+r*cos(p)))) * matrix([-sin(t)*cos(p), cos(t)*cos(p), 0])
    n_p = -(1./(r*r)) * matrix([cos(t)*sin(p), sin(t)*sin(p), -cos(p)])
    return x_t * n_t + x_p * n_p

def KandHwithII(t,p):
    la, v = eig(II(t,p))
    k1, k2 = la[argsort(abs(la))[1:3]]
    return k1*k2, 0.5*(k1+k2)

def IILokal(t,p):
    R = 2.0
    r = 0.5
    x_t = (R+r*cos(p)) * matrix([-sin(t), cos(t), 0.]).T
    x_p = (-r) * matrix([cos(t)*sin(p), sin(t)*sin(p), -cos(p)]).T
    n_t = (1./((R+r*cos(p))*(R+r*cos(p)))) * matrix([-sin(t)*cos(p), cos(t)*cos(p), 0])
    n_p = -(1./(r*r)) * matrix([cos(t)*sin(p), sin(t)*sin(p), -cos(p)])
    S = zeros((2,2))
    S[0,0] = n_t * x_t
    S[0,1] = n_t * x_p
    S[1,0] = n_p * x_t
    S[1,1] = n_p * x_p
    return S

def KandHwithIILokal(t,p):
    la, v = eig(IILokal(t,p))
    k1, k2 = la
    return k1*k2, 0.5*(k1+k2)

def KandH(t,p):
    R = 2.0
    r = 0.5
    return cos(p) / (r * (R+r*cos(p))) , (R+2.*r*cos(p)) / (2.0*r*(R+r*cos(p)))
