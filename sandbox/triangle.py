import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

class triangle(object):

  def __init__(self, p0, p1, p2):
    self._p0 = p0
    self._p1 = p1
    self._p2 = p2
    self._e2 = self._p1 - self._p0
    self._e1 = self._p0 - self._p2
    self._e0 = self._p2 - self._p1

    self.__setCircumcenter()
    self.__setEdgeCircumcenters()

  def __str__(self):
    return "p0: %s\np1: %s\np2: %s\ncircumcenter: %s\n" % (self._p0, self._p1, self._p2, self._circumcenter)

  def __setCircumcenter(self):
    d = 2.0 * norm(cross(self._e1, self._e2))**2
    a1 = - norm(self._e1)**2 * dot(self._e2, self._e0) / d
    a2 = - norm(self._e2)**2 * dot(self._e1, self._e0) / d
    self._circumcenter = self._p0 + a1*self._e2 - a2*self._e1
    print norm(self._circumcenter-self._p0)
    print norm(self._circumcenter-self._p1)
    print norm(self._circumcenter-self._p2)

#  def sinv2(self, a, b):
#    return 2.0 * dot(a, b) * norm(cross(a,b)) / dot(a,a) / dot(b,b)
#
  #def __setCircumcenter(self):
  #  self._circumcenter = self.sinv2(-self._e1, self._e2) * self._p0 + self.sinv2(-self._e2, self._e0) * self._p1 + self.sinv2(-self._e0, self._e1) * self._p2

  #def __setCircumcenter(self):
  #  A = zeros([3,3]);
  #  b = zeros([3,1]);
  #  for i in range(3):
  #    diff01 = self._p0[i] - self._p1[i]
  #    diff02 = self._p0[i] - self._p2[i]
  #    diff12 = self._p1[i] - self._p2[i]
  #    diff01Square = self._p0[i]**2 - self._p1[i]**2
  #    diff02Square = self._p0[i]**2 - self._p2[i]**2
  #    diff12Square = self._p1[i]**2 - self._p2[i]**2
  #    A[0,i] = 2.0 * diff01
  #    A[1,i] = 2.0 * diff02
  #    A[2,i] = 2.0 * diff12
  #    b[0] = b[0] + diff01Square
  #    b[1] = b[1] + diff02Square
  #    b[2] = b[2] + diff12Square
  #  print A
  #  print b
  #  self._circumcenter = linalg.solve(A, b).T[0]
  #  print norm(self._circumcenter-self._p0)
  #  print norm(self._circumcenter-self._p1)
  #  print norm(self._circumcenter-self._p2)

  def __setEdgeCircumcenters(self):
    self._cc1 = self._p2 + 0.5*self._e1
    self._cc0 = self._p1 + 0.5*self._e0
    self._cc2 = self._p0 + 0.5*self._e2

  def plot(self):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = array([self._p0[0], self._p1[0], self._p2[0], self._p0[0]])
    y = array([self._p0[1], self._p1[1], self._p2[1], self._p0[1]])
    z = array([self._p0[2], self._p1[2], self._p2[2], self._p0[2]])
    ax.plot(x, y, z)
    
    ax.plot([self._circumcenter[0], self._p0[0]], [self._circumcenter[1], self._p0[1]], [self._circumcenter[2], self._p0[2]])
    ax.plot([self._circumcenter[0], self._p1[0]], [self._circumcenter[1], self._p1[1]], [self._circumcenter[2], self._p1[2]])
    ax.plot([self._circumcenter[0], self._p2[0]], [self._circumcenter[1], self._p2[1]], [self._circumcenter[2], self._p2[2]])

    xps = [self._circumcenter[0], self._cc0[0], self._cc1[0], self._cc2[0]]
    yps = [self._circumcenter[1], self._cc0[1], self._cc1[1], self._cc2[1]]
    zps = [self._circumcenter[2], self._cc0[2], self._cc1[2], self._cc2[2]]
    ax.plot(xps, yps, zps, "*")
    plt.show();
