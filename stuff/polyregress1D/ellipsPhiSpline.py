from pylab import *

def getNeighborhoodList(n, i, w):
    pos1 = (i-w)%n
    pos2 = (i+w)%n
    if pos1 < pos2:
        return arange(pos1, pos2+1)
    else:
        return append(arange(pos1, n), arange(pos2+1))

def getPMin(ps):
    pMin = zeros(2)
    pMin[0] = min(ps[:,0])
    pMin[1] = min(ps[:,1])
    return pMin

def getPMax(ps):
    pMax = zeros(2)
    pMax[0] = max(ps[:,0])
    pMax[1] = max(ps[:,1])
    return pMax

def orient(p1, p2, p3):
    #print "z= ", p1[0]*(p2[1]-p3[1]) + p2[0]*(p3[1]-p1[1]) + p3[0]*(p1[1]-p2[1])
    return int(sign( p1[0]*(p2[1]-p3[1]) + p2[0]*(p3[1]-p1[1]) + p3[0]*(p1[1]-p2[1]) ))

def intersect(l1, l2):
    if (orient(l1[0], l1[1], l2[0]) ==  orient(l1[0], l1[1], l2[1]) or orient(l2[0], l2[1], l1[0]) ==  orient(l2[0], l2[1], l1[1])):
        return False
    else:
        return True

def isInside(p,V):
    n = len(V)
    pMax = getPMax(V) + pi
    counter = 0
    for i in arange(n-1):
        if intersect([p, pMax], V[i:i+2]):
            #print [p, pMax], V[i:i+2]
            counter += 1 
    if intersect([p, pMax], V[[n-1,0]]):
        counter += 1
        #print [p, pMax], V[[n-1,0]]
    #print counter
    if (counter%2 == 0):
        return False
    else:
        return True


class ellipse:
    def __init__(self, a, b, p0, alpha):
        self.a = a
        self.b = b
        self.p0 = p0
        self.alpha = alpha
    
    def __call__(self, t):
        return array([ self.p0[0] + self.a * cos(t) * cos(self.alpha) - self.b * sin(t) * sin(self.alpha),
                       self.p0[1] + self.a * cos(t) * sin(self.alpha) + self.b * sin(t) * cos(self.alpha)])

    def makeGrid(self, n):
        t = linspace(0, 2.0*pi, n, endpoint=False)
        self.vertices = zeros([n, 2])
        for i in arange(n):
            self.vertices[i] = self(t[i])
        return self.vertices

    def getVertices(self):
        return self.vertices

    def addPlot(self):
        plot(self.vertices[:,0], self.vertices[:,1], "*")

class heart:
    
    def __call__(self, t):
        sint = sin(t)
        return array([ 16.0 * sint * sint * sint,
                       13.0*cos(t) - 5.0*cos(2.0*t) - 2.0*cos(3.0*t) - cos(4.0*t)])

    def makeGrid(self, n):
        t = linspace(0, 2.0*pi, n, endpoint=False)
        self.vertices = zeros([n, 2])
        for i in arange(n):
            self.vertices[i] = self(t[i])
        return self.vertices

    def getVertices(self):
        return self.vertices

    def addPlot(self):
        plot(self.vertices[:,0], self.vertices[:,1], "*")


class elliptSpline:
    def __init__(self, vertices, w):
        self.vertices = array(vertices)
        m = 2 * w +1
        n = len(vertices)
        self.c = zeros([n, 5])
        for i in arange(n):
            locVerts = array(vertices[getNeighborhoodList(n, i, w)])
            if (i > 0):
                self.c[i] = self.c[i-1]
            else:
                #l = locVerts[-1] - locVerts[0]
                #llen = norm(l)
                #lxlen = abs(l[0])
                #self.c[i,0] = llen
                #self.c[i,1] = llen / 2.0
                #self.c[i,2] = arccos(lxlen/llen)
                #self.c[i,3] = locVerts[0,0] + 0.5*l[0]
                #self.c[i,4] = locVerts[0,1] + 0.5*l[1]
                #self.c[i] += 0.1*random(5)
                self.c[i] = [1,2,0,0,0]
            phi = self.getPhiVec(i, locVerts)
            phiAlt = inf * ones(m)
            count = 0
            print "c0   ",  self.c[i]
            print "phi0 ", norm(phi)
            while ( norm(phi-phiAlt) > 1.0e-6 and count < 100 ):
                #print phiAlt, norm(phi-phiAlt), (norm(phi-phiAlt) > 1.0e-2)
                D = self.getParaGradPhiVec(i, locVerts)
                DTD = D.T.dot(D)
                #print i, count, norm(DTD,2)/norm(DTD,-2), norm(phi)
                self.c[i] = self.c[i] - (inv(DTD).dot(D.T).dot(phi))
                phiAlt = phi
                phi = self.getPhiVec(i, locVerts)
                count += 1
                print i, count, norm(DTD,2)/norm(DTD,-2), norm(phi)
                print "c ",  self.c[i]


    def __call__(self, p):
        mini = 0
        minNorm = norm(p - self.vertices[0])
        for i in arange(1, len(self.vertices)):
            normP = norm(p - self.vertices[i])
            if (normP < minNorm):
                minNorm = normP
                mini = i
        val = self.getPhi(mini, p)
        gradPhi = - self.getParaGradPhi(mini, p)[[3,4]]
        e1 = self.vertices[mini] - self.vertices[mini-1]
        e2 = self.vertices[(mini+1)%len(self.vertices)] - self.vertices[mini]
        n1 = -array([e1[1], -e1[0]])
        n2 = -array([e2[1], -e2[0]])
        nAv = 0.5 * (n1 + n2)
        val *= sign(dot(nAv, gradPhi))
        return tanh(5.0*val)
            
    def getPhi(self, i, p):
        a, b, theta, xc, yc = self.c[i]
        xCon = (p[0] - xc) * cos(theta) + (p[1] - yc) * sin(theta)
        yCon = -(p[0] - xc) * sin(theta) + (p[1] - yc) * cos(theta)
        xa = xCon / a
        yb = yCon / b
        return xa * xa + yb * yb - 1.0

    def getPhiVec(self, i, verts):
        n = len(verts)
        rval = zeros(n)
        for k in arange(n):
            rval[k] = self.getPhi(i, verts[k])
        return rval

    def getParaGradPhi(self, i, p):
        a, b, theta, xc, yc = self.c[i]
        xCon = (p[0] - xc) * cos(theta) + (p[1] - yc) * sin(theta)
        yCon = -(p[0] - xc) * sin(theta) + (p[1] - yc) * cos(theta)
        a2 = a * a
        b2 = b * b
        Da = - 2.0 * xCon * xCon / (a2 * a)
        Db = - 2.0 * yCon * yCon / (b2 * b)
        Dtheta = 2.0 * xCon * yCon * ( 1.0/a2 - 1.0/b2 )
        Dxc = 2.0 * (yCon * sin(theta) / b2 - xCon * cos(theta) / a2)
        Dyc = -2.0 * (xCon * sin(theta) / a2 - yCon * cos(theta) / b2)
        return array([Da, Db, Dtheta, Dxc, Dyc])

    def getParaGradPhiVec(self, i, verts):
        n = len(verts)
        rval = zeros([n,5])
        for k in arange(n):
            rval[k,:] = self.getParaGradPhi(i, verts[k])
        return rval

    def addPlot(self, pMin, pMax, n):
        X = linspace(pMin[0], pMax[0], n)
        Y = linspace(pMin[1], pMax[1], n)
        Z = zeros([n,n])
        for  i in arange(len(X)):
            for  j in arange(len(Y)):
                Z[j][i] = self([X[i],Y[j]])
        imshow(Z, interpolation='bilinear', origin='lower', extent=(pMin[0], pMax[0], pMin[1], pMax[1]))
        contour(X,Y,Z,[0.0])


class quadSpline:
    
    def __init__(self, vertices, w):
        self.vertices = array(vertices)
        #meanP = array([mean(vertices[:,0]), mean(vertices[:,1])])
        #meanP = 0.1*random([1,2])
        #print meanP
        m = 2 * w +1
        n = len(vertices)
        self.c = zeros([n, 5])
        conds = zeros(n)

        for i in arange(n):
            locVerts = array(vertices[getNeighborhoodList(n, i, w)])
            xMin = min(locVerts[:,0])
            xMax = max(locVerts[:,0])
            yMin = min(locVerts[:,1])
            yMax = max(locVerts[:,1])
            locVerts[:,0] = (locVerts[:,0] - xMin) / (xMax - xMin)
            locVerts[:,1] = (locVerts[:,1] - yMin) / (yMax - yMin)
            print locVerts
            #locVerts -= mean(locVerts,0)
            #locVerts -= 1.0
            #print locVerts
            #print vertices[getNeighborhoodList(n, i, w)]
            #print locVerts
            A = zeros([m, 5])
            A[:,0:2] = 0.5 * locVerts * locVerts
            A[:,2] = locVerts[:,0] * locVerts[:,1]
            A[:,3:5] = locVerts
            sysMat = A.T.dot(A)
            conds[i] = norm(sysMat,2)/norm(sysMat,-2)
            #print conds[i]
            self.c[i] = inv(sysMat).dot(A.T.dot(ones(m)))
        print "L2Cond:  ", norm(conds)/n
        print "MaxCond: ", max(conds)


    def __call__(self, p):
        mini = 0
        minNorm = norm(p - self.vertices[0])
        for i in arange(1, len(self.vertices)):
            normP = norm(p - self.vertices[i])
            if (normP < minNorm):
                minNorm = normP
                mini = i
        cloc = self.c[mini]
        #print self.vertices[mini]
        #print p
        val = 0.5 * (cloc[0] * p[0] * p[0] + cloc[1] * p[1] * p[1]) + cloc[2] * p[0] * p[1] + cloc[3] * p[0] + cloc[4] * p[1] - 1.0
        
        gradPhi = [cloc[0] * p[0] + cloc[2] * p[1] + cloc[3], cloc[1] * p[1] + cloc[2] * p[0] + cloc[4]]
        e1 = self.vertices[mini] - self.vertices[mini-1]
        e2 = self.vertices[(mini+1)%len(self.vertices)] - self.vertices[mini]
        n1 = -array([e1[1], -e1[0]])
        n2 = -array([e2[1], -e2[0]])
        nAv = 0.5 * (n1 + n2)
        #quiver(self.vertices[mini,0], self.vertices[mini,1], nAv[0], nAv[1])
        #quiver(self.vertices[mini,0], self.vertices[mini,1], n1[0], n1[1])
        #quiver(self.vertices[mini,0], self.vertices[mini,1], n2[0], n2[1])
        #nAv = e1
        val *= sign(dot(nAv, gradPhi))
        #print dot(nAv, gradPhi)
        
        #if ( isInside(p, self.vertices) ):
        #    val *= -1.0

        #if (maxNorm < norm(self.vertices[mini]-self.vertices[maxi])):
        #    val *= -1.0
        #if (abs(val) > 1.0):
        #    val = sign(val)
        return tanh(5.0*val)

    def __call__(self,p):
        phi = 1.0
        for i in arange(len(self.vertices)):
            cloc = self.c[i]
            val = 0.5 * (cloc[0] * p[0] * p[0] + cloc[1] * p[1] * p[1]) + cloc[2] * p[0] * p[1] + cloc[3] * p[0] + cloc[4] * p[1] - 1.0
            if (abs(val) > 1.0):
                phi *= sign(val)
            else:
                phi *= val
        return pow(abs(phi), 1./len(self.vertices))


    def addPlot(self, pMin, pMax, n):
        X = linspace(pMin[0], pMax[0], n)
        Y = linspace(pMin[1], pMax[1], n)
        Z = zeros([n,n])
        for  i in arange(len(X)):
            for  j in arange(len(Y)):
                Z[j][i] = self([X[i],Y[j]])
        imshow(Z, interpolation='bilinear', origin='lower', extent=(pMin[0], pMax[0], pMin[1], pMax[1]))
        contour(X,Y,Z,[0.0])
        


#el = ellipse(1,2, [0,0], pi/16.0)
#el = ellipse(1,2, [0,0], 0.0)
el = heart()
V = el.makeGrid(100)
el.addPlot()
#for i in arange(len(V)-1):
#    plot(V[i:i+2, 0], V[i:i+2, 1], "->")

spline = quadSpline(V, 3)
#spline.addPlot(getPMin(V), getPMax(V), 20)

#spline = elliptSpline(V, 3)
#spline.addPlot(getPMin(V), getPMax(V), 20)

show()
