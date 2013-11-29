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

    xps = [self._circumcenter[0], self._cc0[0], self._cc1[0], self._cc2[0]]
    yps = [self._circumcenter[1], self._cc0[1], self._cc1[1], self._cc2[1]]
    zps = [self._circumcenter[2], self._cc0[2], self._cc1[2], self._cc2[2]]
    ax.plot(xps, yps, zps, "*")
    plt.show();
