from pylab import *

dt = float32

#
#         
#         w_3 w_2
#          \   |      W_1
#           \  |    _/
#            \ |  _/  
#             \| /
#   w_4 ------ v ------- w_0
#           .   \
#            ... \_
#                  \
#                   w_n
#
#
#
#

#def csd_vols(v1, v2, wL) :
#    edge = v2 - v1
#    L1_edge = wL - v1
#    L2_edge = wL - v2
#
#    TL_vol = 0.5 * abs(edge[0] * L1_edge[1] - edge[1] * L1_edge[0])
#
#    cL = v1 + (dot(edge, edge) * dot(L1_edge, L2_edge) * L1_edge - dot(L1_edge, L1_edge) * dot(edge, L2_edge) * edge) / (8.0 * TL_vol * TL_vol)
#
#    hh12 = norm(cL - 0.5 * (v2 + v1))
#    hh1L = norm(cL - 0.5 * (wL + v1))
#
#    vvol1 = 0.25 * (norm(edge) * hh12 + norm(L1_edge) * hh1L)
#
#    return vvol1, hh12, hh1L # dual vol at v1 , dual len at (v2 - v1), dual len at (wL - v1)

def csd_vols(v, w1, w2) :
    e1 = w1 - v
    e2 = w2 -v
    e21 = w2 - w1
    
    h1 = norm(e1)
    h2 = norm(e2)
    h21 = norm(e21)

    cg21 = dot(e1, e2) / (h1 * h2)
    cg1  = dot(e2, e21) / (h2 * h21)

    hh1 = h1 * ( h2 / h1 - cg21) / sqrt( 1.0 - cg21**2)
    hh2 = h2 * ( h21 / h2 - cg1) / sqrt( 1.0 - cg1**2)

    vvol = 0.25 * (h1 * hh1 + h2 * hh2)

    return vvol, hh1, hh2 # dual vol at v , dual len at (w1 - v), dual len at (w2 - v)

def common_quantities(v, ws):
    N_exvertices = len(ws)
    hv_vol = dt(0.0)
    hes_len = zeros(N_exvertices, dtype=dt)
    for i in arange(N_exvertices):
        vvol, hh1, hh2 = csd_vols(v, ws[i], ws[(i+1)%N_exvertices])
        hes_len[i] += hh1
        hes_len[(i+1)%N_exvertices] += hh2
        hv_vol += vvol
    return hv_vol, hes_len # voronoi vol |*v| , voronoi lengths |*e|
    
    
class DivMesh:
        
    def __init__(self, v, hv_vol, ws, hes_len):
        self.v = array(v, dtype=dt)  # center vertex
        self.hv_vol = dt(hv_vol) # associated dual area resp. to the center point (common DEC: voronoi area)
        self.ws = array(ws, dtype=dt) # exterior vertices
        self.hes_len = array(hes_len, dtype=dt) # associated dual lenght resp. to the edges (w_i - v) (common DEC: |*e_i|)

        self.h = norm(self.ws[0] - self.v)
    
    @classmethod
    def create_equidistance_mesh_common(cls, N_exvertices, h=1, v=[0,0]):
        h = dt(h)
        v = array(v, dtype=dt)

        angle_step = 2. * dt(pi) / N_exvertices #FIXME: pi of dt type?
        ws = zeros((N_exvertices,2), dtype=dt)
        for i in arange(N_exvertices):
            angle = i * angle_step
            ws[i] = v + h * array([cos(angle), sin(angle)], dtype=dt)
        
        hv_vol, hes_len = common_quantities(v, ws) 

        return cls(v, hv_vol, ws, hes_len)

    @classmethod
    def create_disturbed_equidistance_mesh_common(cls, N_exvertices, h=1, v=[0,0], seed=42, disturb_fac=1.0, disturb_kind='uniform', disturb_impact='all', disturb_direction='both'):
        # disturb_kind='uniform'|'static'
        # disturb_impact='all'|'single'
        # disturb_direction='angle'|'radius'|'both'
        np.random.seed(seed)
        h = dt(h)
        v = array(v, dtype=dt)
        
        if (disturb_kind == 'uniform'):
            rn = lambda : disturb_fac * (2.*rand() - 1.) # in (-fac,fac)
        elif (disturb_kind == 'static'):
            rn = lambda: disturb_fac
        else:
            raise ValueError("'disturb_kind' must be 'uniform' or 'static'")

        if (disturb_direction=='angle'):
            rna = lambda : rn()
            rnr = lambda : 0.
        elif (disturb_direction=='radius'):
            rna = lambda : 0.
            rnr = lambda : rn()
        elif (disturb_direction=='both'):
            rna = lambda : rn()
            rnr = lambda : rn()
        else:
            raise ValueError("'disturb_direction' must be 'angle', 'radius' or 'both'")

        if (not(disturb_impact=='all' or disturb_impact=='single')):
            raise ValueError("'disturb_impact' must be 'single' or 'all'")
    
        angle_step = 2. * dt(pi) / N_exvertices #FIXME: pi of dt type?
        ws = zeros((N_exvertices,2), dtype=dt)
        for i in arange(N_exvertices):
            angle = (i + rna()) * angle_step
            ws[i] = v + (h + rnr()) * array([cos(angle), sin(angle)], dtype=dt)
            if (disturb_impact=='single'):
                rn = lambda : 0.0
        
        hv_vol, hes_len = common_quantities(v, ws) 

        return cls(v, hv_vol, ws, hes_len)


    @classmethod
    def create_scaled_mesh_by_factor(cls, divMesh, scale_fac):
        #length dimensions
        v = scale_fac * divMesh.v
        ws = scale_fac * divMesh.ws
        hes_len = scale_fac * divMesh.hes_len
        #area dimensions
        hv_vol = scale_fac * scale_fac * divMesh.hv_vol
        return cls(v, hv_vol, ws, hes_len)

    def make_consistence(self, order):
        self.hv_vol = dt(pi) * self.h**2 / dt(4)

        if (order == 1):
            wm = self.W1_mat()
            rv = self.R1_vec()
        elif (order ==2):
            wm = self.W2_mat()
            rv = self.R2_vec()
        else:
            raise ValueError("'order' must be in {1,2}")
            
        wmt = transpose(wm)
        #cv = solve( dot(wmt,wm) , dot(wmt, rv)) #NR
        cv = dot(wmt , solve( dot(wm,wmt) , rv)) #NE
        for i in arange(len(self.hes_len)):
            self.hes_len[i] = cv[i] * norm(self.v - self.ws[i])


    def divergence(self, form):
        divvol = dt(0.0)
        N= len(self.ws)
        for i in arange(N):
            edge = self.ws[i] - self.v
            divvol += (self.hes_len[i] / norm(edge)) * form.integrate(edge, self.v)
        return divvol / self.hv_vol

    def W1_mat(self):
        N= len(self.ws)
        wmat = zeros((5,N), dtype=dt)
        for i in arange(N):
            wmat[0][i] = self.ws[i][0]
            wmat[1][i] = self.ws[i][1]
            wmat[2][i] = self.ws[i][0]**2
            wmat[3][i] = self.ws[i][1]**2
            wmat[4][i] = self.ws[i][0] * self.ws[i][1]
        return wmat

    def W2_mat(self):
        N= len(self.ws)
        wmat = zeros((9,N), dtype=dt)
        for i in arange(N):
            wmat[0][i] = self.ws[i][0]
            wmat[1][i] = self.ws[i][1]
            wmat[2][i] = self.ws[i][0]**2
            wmat[3][i] = self.ws[i][1]**2
            wmat[4][i] = self.ws[i][0] * self.ws[i][1]
            wmat[5][i] = self.ws[i][0]**3
            wmat[6][i] = self.ws[i][1]**3
            wmat[7][i] = (self.ws[i][0]**2) * self.ws[i][1]
            wmat[8][i] = self.ws[i][0] * (self.ws[i][1]**2)
        return wmat

    def R1_vec(self):
        return 2.0*self.hv_vol*array([0,0,1,1,0],dtype=dt)

    def R2_vec(self):
        return 2.0*self.hv_vol*array([0,0,1,1,0, 0,0,0,0],dtype=dt)

    def C_vec(self):
        N= len(self.ws)
        cvec = zeros(N, dtype=dt)
        for i in arange(N):
            cvec[i] = self.hes_len[i] / norm(self.v - self.ws[i])
        return cvec

    def plot_mesh(self):
        edge_color = 'b'
        N= len(self.ws)
        for i in arange(N):
            plot([self.v[0],self.ws[i][0]], [self.v[1],self.ws[i][1]], color=edge_color)
            plot([self.ws[i-1][0],self.ws[i][0]], [self.ws[i-1][1],self.ws[i][1]], color=edge_color)


# u := [ax, ay] * (x+1)^2 * (y+1)^2  , ax,ay in R
class ExampleForm_1:
    
    def __init__(self, ax, ay):
        self.ax = dt(ax)
        self.ay = dt(ay)

        self.xyfun = lambda x, y: (x+1.0) * (x+1.0) * (y+1.0) * (y+1.0) * array([self.ax,self.ay],dtype=dt)
        self.xyfun_div = lambda x, y: 2.0 * (x+1.0) * (y+1.0) * (self.ax * (y+1.0) + self.ay * (x+1.0))
        #first order error derivative for Nv=3
        #self.xyfun_DNv3 = lambda x, y: (4.0 * self.ay * (x+1.0) * (y+1.0) + self.ax * ( x - y ) * ( 2. + x + y)) / 3.0
        
    def integrate(self, edge, v1):
        v11p1 = v1[0] + 1.0
        v12p1 = v1[1] + 1.0
        return (1.0 / 30.) * (self.ax * edge[0] + self.ay * edge[1]) * (10.0 * v11p1 * v11p1 * (edge[1]*edge[1] + 3.0 * edge[1] * v12p1 + 3.0 * v12p1 * v12p1) + 5.0 * edge[0] * v11p1 * (3.0 * edge[1]*edge[1] + 8.0 * edge[1] * v12p1 + 6.0 * v12p1 * v12p1) + edge[0] * edge[0] * (6.0 * edge[1]*edge[1] + 15.0 * edge[1] * v12p1 + 10.0 * v12p1 * v12p1))

# u := [ax, ay] = const
class ExampleForm_2:
    
    def __init__(self, ax, ay):
        self.ax = dt(ax)
        self.ay = dt(ay)

        self.xyfun = lambda x, y: array([self.ax,self.ay],dtype=dt)
        self.xyfun_div = lambda x, y: 0
        
    def integrate(self, edge, v1):
        return edge[0]*self.ax + edge[1]*self.ay



class Test_Div:
    
    def __init__(self, meshes, form, n_refinements):
        #self.DNv3AtV = form.xyfun_DNv3(meshes[0].v[0], meshes[0].v[1])
        self.nm = len(meshes)
        self.nr = n_refinements
        self.hs = []
        self.errs = []
        self.res0s_rel1 = []
        self.res0s_relvvol = []
        self.res1s_rel2 = []
        self.res1s_relvvol = []
        for mesh in meshes:
            self.hs.append([])
            self.errs.append([])
            self.res0s_rel1.append([])
            self.res0s_relvvol.append([])
            self.res1s_rel2.append([])
            self.res1s_relvvol.append([])
            for scale in concatenate(([1.],ones(self.nr)*0.5)):
                mesh = DivMesh.create_scaled_mesh_by_factor(mesh, scale)
                self.hs[-1].append(mesh.h)

                print form.xyfun(mesh.v[0], mesh.v[1])
                div = form.xyfun_div(mesh.v[0], mesh.v[1])
                divh = mesh.divergence(form)
                self.errs[-1].append(abs(div - divh))
                
                wc2mr = dot(mesh.W2_mat(),mesh.C_vec()) - mesh.R2_vec()
                #self.res1s[-1].append(norm(wc1mr))
                self.res0s_rel1[-1].append(norm(wc2mr[0:2])/(mesh.h**1))
                self.res0s_relvvol[-1].append(norm(wc2mr[0:2])/(mesh.hv_vol))
                self.res1s_rel2[-1].append(norm(wc2mr[2:5])/(mesh.h**2))
                self.res1s_relvvol[-1].append(norm(wc2mr[2:5])/(mesh.hv_vol))

    def print_results(self):
        for im in arange(self.nm):
            print ""
            for ir in arange(self.nr+1):
                print "h = ", self.hs[im][ir], "; Error = ",  self.errs[im][ir], "; Res0_Rel1 = ",  self.res0s_rel1[im][ir], "; Res0_RelVVol = ",  self.res0s_relvvol[im][ir], "; Res1_Rel2 = ",  self.res1s_rel2[im][ir], "; Res1_RelVVol = ",  self.res1s_relvvol[im][ir]
                #print self.errs[im][ir] / self.hs[im][ir], " ; ",  self.DNv3AtV, " ; ", abs(self.errs[im][ir] / self.hs[im][ir] - self.DNv3AtV)




# main
#mesh = DivMesh.create_equidistance_mesh_common(5)
#print mesh.hes_len
#mesh = DivMesh.create_disturbed_equidistance_mesh_common(5, seed=42, disturb_fac=2.24, disturb_kind="static", disturb_impact='single', disturb_direction='radius') #O(h^-1)
#mesh = DivMesh.create_disturbed_equidistance_mesh_common(5, seed=42, disturb_fac=2.23, disturb_kind="static", disturb_impact='single', disturb_direction='radius') #O(h)
#mesh = DivMesh.create_disturbed_equidistance_mesh_common(5, seed=42, disturb_fac=0.1, disturb_kind="static", disturb_impact='single', disturb_direction='radius') #O(h)
#mesh = DivMesh.create_disturbed_equidistance_mesh_common(5, seed=42, disturb_fac=0.4, disturb_kind="uniform", disturb_impact='all', disturb_direction='both')
#mesh = DivMesh.create_disturbed_equidistance_mesh_common(3, seed=42, disturb_fac=0.5, disturb_kind="static", disturb_impact='all', disturb_direction='angle')

fac = 0.0
mesh = DivMesh.create_disturbed_equidistance_mesh_common(5, seed=12, disturb_fac=fac, disturb_kind="uniform", disturb_impact='all', disturb_direction='both')
#mesh = DivMesh.create_scaled_mesh_by_factor(mesh, 100000)
#mesh = DivMesh.create_disturbed_equidistance_mesh_common(7, seed=42, disturb_fac=fac, disturb_kind="static", disturb_impact='single', disturb_direction='radius')


#mesh.make_consistence(1)
print mesh.h
print mesh.hes_len

mesh.plot_mesh()
form = ExampleForm_2(1,-4.2)

meshes = [mesh]
for i in arange(0):
    fac *= 0.5
    meshes.append( DivMesh.create_disturbed_equidistance_mesh_common(7, seed=42, disturb_fac=fac, disturb_kind="uniform", disturb_impact='all', disturb_direction='both'))
    meshes[-1].plot_mesh()
test = Test_Div(meshes, form, 30)
test.print_results()

plt.axes().set_aspect('equal')
show()
