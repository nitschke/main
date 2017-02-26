#import matplotlib.pyplot as plt
from pylab import *
dt = float128
#           wL 
#          ^ ^
#      L1 / ^ \ L2  TL
#        /  |  \
#      v1------>v2
#        \ e|  /
#      R1 \ | / R2  TR
#          v v 
#           wR
#
class RotRotMesh:
    
    # init "unweighted" (common circumcentered DEC)
    def __init__(self, v1, v2, wL, wR):
        self.v1 = array(v1, dtype=dt)
        self.v2 = array(v2, dtype=dt)
        self.wL = array(wL, dtype=dt)
        self.wR = array(wR, dtype=dt)

        self.edge = self.v2 - self.v1
        self.L1_edge = self.wL - self.v1
        self.L2_edge = self.wL - self.v2
        self.R1_edge = self.wR - self.v1
        self.R2_edge = self.wR - self.v2
        
        self.TL_vol = 0.5 * abs(self.edge[0] * self.L1_edge[1] - self.edge[1] * self.L1_edge[0])
        self.TR_vol = 0.5 * abs(self.edge[0] * self.R1_edge[1] - self.edge[1] * self.R1_edge[0])
        
        #circumcenters: see diploma (1.29)
        self.cL = self.v1 + (dot(self.edge, self.edge) * dot(self.L1_edge, self.L2_edge) * self.L1_edge - dot(self.L1_edge, self.L1_edge) * dot(self.edge, self.L2_edge) * self.edge) / (8.0 * self.TL_vol * self.TL_vol)
        self.cR = self.v1 + (dot(self.edge, self.edge) * dot(self.R1_edge, self.R2_edge) * self.R1_edge - dot(self.R1_edge, self.R1_edge) * dot(self.edge, self.R2_edge) * self.edge) / (8.0 * self.TR_vol * self.TR_vol)
        
        self.c = self.v1 + 0.5 * self.edge
        
        self.h = norm(self.edge)

    # init "unweighted" (common circumcentered DEC) with equilateral triangles
    @classmethod
    def create_equilateral_triangles(cls, v1, v2):
        edge =  array(v2, dtype=dt) - array(v1, dtype=dt)
        edge_rotated = array([-edge[1], edge[0]], dtype=dt)
        fac = sqrt(3.) / 2. # he / h
        return cls(v1, v2, v1 + 0.5 * edge + fac * edge_rotated, v1 + 0.5 * edge - fac * edge_rotated)

    # init "unweighted" (common circumcentered DEC) with flipped triangles, then displace wR by wR_offset
    @classmethod
    def create_flipped_triangles(cls, v1, v2, wL, wR_offset=array([0.,0.], dtype=dt)):
        v1 = array(v1, dtype=dt)
        v2 = array(v2, dtype=dt)
        wL = array(wL, dtype=dt)
        wR_offset = array(wR_offset, dtype=dt)
        edge = v2 - v1
        h = norm(edge)
        flipmap = (2. / (h * h)) * outer(edge, edge) - eye(2)
        wR = v1 + dot(flipmap, wL - v1) + wR_offset
        return cls(v1, v2, wL, wR)
    
    def optimize_duallength(self):
        dualedge = self.cL - self.cR
        h_dualedge = norm(dualedge)
        scalefac = abs(dot(self.wL - self.wR, dualedge)) / (3. * h_dualedge * h_dualedge)
        dualedge_L = self.cL - self.c
        dualedge_R = self.cR - self.c
        self.cL = self.c + scalefac * dualedge_L
        self.cR = self.c + scalefac * dualedge_R

    def scale(self, scale_fac):
        self.v1 *= scale_fac 
        self.v2 *= scale_fac
        self.wL *= scale_fac
        self.wR *= scale_fac
        self.edge *= scale_fac
        self.L1_edge *= scale_fac
        self.L2_edge *= scale_fac
        self.R1_edge *= scale_fac
        self.R2_edge *= scale_fac
        self.cL *= scale_fac
        self.cR *= scale_fac
        self.h *= scale_fac
        self.c *= scale_fac


        self.TL_vol *= scale_fac * scale_fac
        self.TR_vol *= scale_fac * scale_fac


    def plot(self):
        vertices_line = transpose([self.v1,self.wL,self.v2,self.wR,self.v1,self.v2])
        plot(vertices_line[0],vertices_line[1], "-o")
        annotate("$v_1$", xy=self.v1, xytext=(self.v1 + 0.1 * self.h * array([-1.,0.])))
        annotate("$v_2$", xy=self.v2, xytext=(self.v2 + 0.1 * self.h * array([1.,0.])))
        annotate("$v_L$", xy=self.wL, xytext=(self.wL + 0.1 * self.h * array([0.,1.])))
        annotate("$v_R$", xy=self.wR, xytext=(self.wR + 0.1 * self.h * array([0.,-1.])))

        dual_line = transpose([self.cR,self.cL])
        plot(dual_line[0],dual_line[1], "-o")
        annotate("$\star T_L$", xy=self.cL, xytext=(self.cL + 0.1 * self.h * array([0.,1.])))
        annotate("$\star T_R$", xy=self.cR, xytext=(self.cR + 0.1 * self.h * array([0.,-1.])))

    def print_pgf_coords(self):
        print "\coordinate (V1) at (", self.v1[0], ",", self.v1[1], ");"
        print "\coordinate (V2) at (", self.v2[0], ",", self.v2[1], ");"
        print "\coordinate (VL) at (", self.wL[0], ",", self.wL[1], ");"
        print "\coordinate (VR) at (", self.wR[0], ",", self.wR[1], ");"
        print "\coordinate (CL) at (", self.cL[0], ",", self.cL[1], ");"
        print "\coordinate (CR) at (", self.cR[0], ",", self.cR[1], ");"
        print "\coordinate (C) at (", self.c[0], ",", self.c[1], ");"


# u := [ax, ay] * (x+1)^2 * (y+1)^2  , ax,ay in R
class ExampleForm_1:
    
    def __init__(self, ax, ay):
        self.xyfun = lambda x, y: (x+1.0) * (x+1.0) * (y+1.0) * (y+1.0) * array([ax,ay])
        self.xyfun_rotrot = lambda x, y: array([2.0 * ax * (x+1.0) * (x+1.0) - 4.0 * ay * (x+1.0) * (y+1.0), 2.0 * ay * (y+1.0) * (y+1.0) - 4.0 * ax * (x+1.0) * (y+1.0)])
        
        self.ax = ax
        self.ay = ay

    def integrate(self, edge, v1):
        v11p1 = v1[0] + 1.0
        v12p1 = v1[1] + 1.0
        return (1.0 / 30.) * (self.ax * edge[0] + self.ay * edge[1]) * (10.0 * v11p1 * v11p1 * (edge[1]*edge[1] + 3.0 * edge[1] * v12p1 + 3.0 * v12p1 * v12p1) + 5.0 * edge[0] * v11p1 * (3.0 * edge[1]*edge[1] + 8.0 * edge[1] * v12p1 + 6.0 * v12p1 * v12p1) + edge[0] * edge[0] * (6.0 * edge[1]*edge[1] + 15.0 * edge[1] * v12p1 + 10.0 * v12p1 * v12p1))

    def integrate_rotrot(self, edge, v1):
        v12 = v1[1]
        v11 = v1[0]
        e1 = edge[0]
        e2 = edge[1]
        t1 = 1. + v12
        t2 = 1. + v11
        t4 = 3. * t2 * (2. + e2 + 2. * v12)
        t5 = t4 + e1 * (3. + 2. * e2 + 3. * v12)
        t6 = e2 * e2 + 3. * e2 * t1 + 3. * t1 * t1
        t7 = e1 * e1 + 3. * e1 * t2 + 3. * t2 * t2
        return (2./3.) * (e2 * (-self.ax * t5 + self.ay * t6) + e1 * (-self.ay * t5 + self.ax * t7))

class Rotrot_Laplace:

    def __init__(self, mesh, form, duallength_optimization=False):
        hh = norm(mesh.cL - mesh.cR)
        if (duallength_optimization):
            hh *= (mesh.wL[1] - mesh.wR[1])/ ( 3. * norm(mesh.cL - mesh.cR))
            #hh *= abs(dot(mesh.wL - mesh.wR, mesh.cL - mesh.cR)) / (3. * hh * hh)
        form_e = form.integrate(mesh.edge, mesh.v1)
        form_eL1 = form.integrate(mesh.L1_edge, mesh.v1)
        form_eR1 = form.integrate(mesh.R1_edge, mesh.v1)
        form_eL2 = form.integrate(mesh.L2_edge, mesh.v2)
        form_eR2 = form.integrate(mesh.R2_edge, mesh.v2)

        form_e_rotrot = form.integrate_rotrot(mesh.edge, mesh.v1)
        
        self.rotrot_e = (mesh.h/hh) * ((form_eR1 - form_eR2 - form_e) / mesh.TR_vol + (form_eL1 - form_eL2 - form_e) / mesh.TL_vol)
        self.error = form_e_rotrot/ mesh.h - self.rotrot_e / mesh.h

        trptl = mesh.TR_vol + mesh.TL_vol
        form_e_sol = - (mesh.TR_vol * mesh.TL_vol * hh / (trptl * mesh.h)) * form_e_rotrot + (mesh.TL_vol / trptl) * (form_eR1 - form_eR2) + (mesh.TR_vol / trptl) * (form_eL1 - form_eL2)

        self.error_sol = form_e/mesh.h - form_e_sol/mesh.h


class Test_Rotrot:

    def __init__(self, meshes, example, n_refinements):
        self.nm = len(meshes)
        self.nr = n_refinements
        self.erry_prop = []
        self.errx_prop = []
        self.hs = []
        self.errInts = []
        self.errForms = []
        self.errSolInts = []
        for mesh in meshes:
            dualedge = mesh.cL - mesh.cR
            hh = norm(dualedge)
            unitvecP = mesh.edge / mesh.h
            unitvecD = dualedge / hh
            sol_c = example.xyfun_rotrot(mesh.c[0], mesh.c[1])
            self.erry_prop.append((mesh.wL[0] - mesh.wR[0]) * dot(sol_c, unitvecD) / ( 3. * hh))
            self.errx_prop.append((1. - (mesh.wL[1] - mesh.wR[1])/ ( 3. * hh)) * dot(sol_c, unitvecP))
            self.hs.append([])
            self.errInts.append([])
            self.errForms.append([])
            self.errSolInts.append([])
            for scale in concatenate(([1.],ones(self.nr)*0.5)):
                mesh.scale(scale)
                decOp = Rotrot_Laplace(mesh, example, duallength_optimization=False)
                self.hs[-1].append(mesh.h)
                self.errInts[-1].append(decOp.error)
                self.errForms[-1].append(dot(sol_c, unitvecP) - decOp.rotrot_e / mesh.h)
                self.errSolInts[-1].append(decOp.error_sol)

    def print_results(self):
        for im in arange(self.nm):
            print "\nPropagated Errors: (", self.erry_prop[im], " , " , self.errx_prop[im], "); Sum: ", self.erry_prop[im] + self.errx_prop[im]
            for ir in arange(self.nr+1):
                print "h = ", self.hs[im][ir], "; ErrorInt = ",  self.errInts[im][ir], "; ErrorForm = ",  self.errForms[im][ir], "; ErrorSolInt = ",  self.errSolInts[im][ir]

    def plot_results(self, names):
        f, axarr = plt.subplots(2, sharex=True)
        for im in arange(self.nm):
            errsum = self.erry_prop[im] + self.errx_prop[im]
            axarr[0].semilogx(self.hs[im], self.errInts[im], label=names[im], marker="*")
            axarr[0].axhline(errsum, color="k", linestyle="--");
            axarr[1].loglog(self.hs[im], abs(self.errInts[im] - errsum), label=names[im], marker="*")
            #axarr[1].loglog(self.hs[im], abs(array(self.errSolInts[im])), label=names[im]+" (SE)", marker="*")
        axarr[0].legend(loc="upper left", ncol = self.nm)
        axarr[1].legend(loc="upper left", ncol = self.nm)
        axarr[1].set_xlabel("h")
        axarr[0].set_ylabel("$\operatorname{Err}_{h}^{\Delta^\operatorname{Rr}}(e)$")
        axarr[1].set_ylabel("$\operatorname{Err}_{h}^{\Delta^\operatorname{Rr}}(e) - R_{0}^{\Delta^\operatorname{Rr}}\left( c(e) \\right)$")

# main
#mesh = RotRotMesh([-1,0],[1,0],[0.5,1.5],[0,-2])
#mesh = RotRotMesh.create_flipped_triangles([-1,0],[1,0],[0.5,1.5], [0,0.5])
#mesh.optimize_duallength()
#mesh = RotRotMesh.create_equilateral_triangles([-1,0],[1,0])
#mesh.print_pgf_coords()

#mesh.plot()
#n_refinements = 10
#example = ExampleForm_1(1.0, 1.0)
#errx_prop = (mesh.wL[0] - mesh.wR[0]) * example.xyfun_rotrot(0.0,0.0)[1] / ( 3. * norm(mesh.cL - mesh.cR))
#erry_prop = (1. - (mesh.wL[1] - mesh.wR[1])/ ( 3. * norm(mesh.cL - mesh.cR))) * example.xyfun_rotrot(0.0,0.0)[0] 
#print "propagated x-error" , errx_prop
#print "propagated y-error" , erry_prop
#print "propagated error" , errx_prop + erry_prop
#for scale in concatenate(([1.],ones(n_refinements)*0.5)):
#    mesh.scale(scale)
#    #mesh.plot()
#    decOp = Rotrot_Laplace(mesh, example, duallength_optimization=False)
#    errInt = decOp.error
#    errForm = -(decOp.rotrot_e / mesh.h - example.xyfun_rotrot(0.0,0.0)[0])
#    print "h = ", mesh.h, "; ErrorInt = ",  errInt, "; ErrorForm= ",  errForm
#
#    ## test integration
#    #errInt = example.integrate(mesh.edge,mesh.v1) / mesh.h - example.xyfun(0.0,0.0)[0] #dot(example.xyfun(0.0,0.0),mesh.edge)/mesh.h
#    #errInt_rotrot = example.integrate_rotrot(mesh.edge,mesh.v1) / mesh.h - example.xyfun_rotrot(0.0,0.0)[0] #dot(example.xyfun_rotrot(0.0,0.0),mesh.edge)/mesh.h
#    #print "h = ", mesh.h, "; ErrorInt = ",  errInt, "; ErrorInt_Rotrot = ",  errInt_rotrot

example = ExampleForm_1(1.0, 1.0)
K1 = RotRotMesh([-1,0],[1,0],[0.5,1.5],[0,-2])
K2 = RotRotMesh.create_flipped_triangles([-1,0],[1,0],[0.5,1.5], [0,0.5])
K3 = RotRotMesh([-1,0],[1,0],[0.5,1.5],[0,-2])
K3.optimize_duallength()
K4 = RotRotMesh.create_flipped_triangles([-1,0],[1,0],[0.5,1.5], [0,0.5])
K4.optimize_duallength()
K5 = RotRotMesh.create_equilateral_triangles([-1,0],[1,0])
test = Test_Rotrot([K1,K2,K3,K4, K5], example, 21)
test.print_results()
test.plot_results(["$\mathcal{K}_1$", "$\mathcal{K}_2$", "$\mathcal{K}_3$", "$\mathcal{K}_4$","$\mathcal{K}_5$" ])

#plt.axes().set_aspect('equal')
show()

