from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from numpy import *


class Dual_triangulation:
    def __init__(self, tri):
        self.tri = tri

        self.n_vertices = len(tri.points)
        self.n_faces = len(tri.simplices)

        self.w = ones(self.n_vertices)
        
        self.c_faces = zeros([self.n_faces,2])
        for i in arange(self.n_faces):
            face = tri.simplices[i]
            e01 = tri.points[face][1] - tri.points[face][0]
            e12 = tri.points[face][2] - tri.points[face][1]
            e20 = tri.points[face][0] - tri.points[face][2]
            l01 = linalg.norm(e01)
            l12 = linalg.norm(e12)
            l20 = linalg.norm(e20)
            s = 0.5 * (l01 + l12 + l20)
            area = sqrt(s * (s - l01) * (s - l12) * (s - l20))
            self.c_faces[i] = tri.points[face][0] + ((l01*l01 + self.w[face[0]] - self.w[face[1]])*cross([0,0,1],e20)[0:2] + (l20*l20 + self.w[face[0]] - self.w[face[2]])*cross([0,0,1],e01)[0:2]) / area / 4.0




random.seed(42)
nrow= 5
h=1./(nrow-1.0)
n= nrow*nrow
points= zeros([n,2])
for i in arange(nrow):
    for j in arange(nrow):
        points[i+nrow*j,:] = [i*h,j*h]
        if (i != 0 and i != nrow-1 and j != 0 and j != nrow-1):
             points[i+nrow*j,:] += 0.1*(2*random.rand(2) - 1)


plt.plot(points[:,0], points[:,1], 'o')


tri = Delaunay(points)
plt.triplot(points[:,0], points[:,1], tri.simplices.copy())


dualtri = Dual_triangulation(tri)
plt.plot(dualtri.c_faces[:,0], dualtri.c_faces[:,1], 'o')

plt.axes().set_aspect('equal')
plt.show()
