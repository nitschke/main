from pylab import *

def II(u,v):
    x_u = matrix([cos(u)*cos(v), cos(u)*sin(v), -sin(u)])
    x_v = (1./sin(u)) * matrix([-sin(v), cos(v), 0])
    return x_u.T * x_u + x_v.T * x_v

def II2(u,v):
    x_u = matrix([cos(u)*cos(v), cos(u)*sin(v), -sin(u)])
    xov = sin(u)*matrix([-sin(v), cos(v), 0])
    return x_u.T * x_u + xov.T * xov

def II3(u,v):
    x_u = matrix([cos(u)*cos(v), cos(u)*sin(v), -sin(u)])
    x_v = (1./sin(u)) * matrix([-sin(v), cos(v), 0])
    xov = sin(u)*matrix([-sin(v), cos(v), 0])
    return x_u.T * x_u + xov.T * x_v

