#include "Dec.h"
#include "io/VtkVectorWriter.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

// 1-form -> Rot(z)
class Alpha : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Alpha() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    double con2 = 1.0; //contra coords
    WorldVector<double> conBasis2; //contra basis vectors
    conBasis2[0] = -x[1];
    conBasis2[1] = x[0];
    conBasis2[2] = 0.0;
    return (con2 * conBasis2) * vec;
  }
};

// <alpha,[p,q]>
class Alpha_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Alpha_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return 2.*(p2*q1 - 1.*p1*q2)*atan((-1. + p1*q1 + p2*q2 + p3*q3)*pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + 
       pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5))*
   pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + 
     pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5);
  }
};


// Laplace-Beltrami of alpha -> -2 * alpha
class LbAlpha : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LbAlpha() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    double con2 = 1.0; //contra coords
    WorldVector<double> conBasis2; //contra basis vectors
    conBasis2[0] = -x[1];
    conBasis2[1] = x[0];
    conBasis2[2] = 0.0;
    return -2.0 * (con2 * conBasis2) * vec;
  }
};

// <lbalpha,[p,q]>
class LbAlpha_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LbAlpha_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return -4.*(p2*q1 - 1.*p1*q2)*atan((-1. + p1*q1 + p2*q2 + p3*q3)*pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + 
       pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5))*
   pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + 
     pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5);
  }
};


// exact 1-form : gamma = d(xyz)
class DXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    WorldVector<double> gradF;
    for (int i = 0; i < 3; i++) {
      int ii = (i+1)%3;
      int iii = (i+2)%3;
      gradF[i] = x[ii] * x[iii] * (1.0 - 3.0*x[i]*x[i]);
    }
    return  gradF * vec;
  }
};

// <dxyz,[p,q]>
class DXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
  }
};

// <Lcbdxyz,[p,q]> -> -12 * dxyz
class LcbDXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LcbDXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  -12.*(q[0]*q[1]*q[2] - p[0]*p[1]*p[2]);
  }
};

// rot(x*y*z)
class RotXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    WorldVector<double> rot;
    for (int i = 0; i < 3; i++) {
      int ii = (i+1)%3;
      int iii = (i+2)%3;
      rot[i] = x[i] * (x[ii]*x[ii] - x[iii]*x[iii]);
    }
    return  rot * vec;
  }
};



// Laplace-Beltrami of rot(x*y*z) -> -12*rot(x*y*z)
class LbRotXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  LbRotXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    WorldVector<double> rot;
    for (int i = 0; i < 3; i++) {
      int ii = (i+1)%3;
      int iii = (i+2)%3;
      rot[i] = x[i] * (x[ii]*x[ii] - x[iii]*x[iii]);
    }
    return  -12.0 * (rot * vec);
  }
};




// e.g. Laplace-Beltrami of exact forms or Laplace-CoBeltrami of co-exact forms
class Zero : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Zero() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    return 0.0;
  }
};


//projection
class Proj : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

  public:
  Proj() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    return x / sqrt(x * x);
  }
};

//jacobi matrix of projection
class JProj : public AbstractFunction<WorldMatrix<double>, WorldVector<double> > {

  public:
  JProj() : AbstractFunction<WorldMatrix<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldMatrix<double> operator()(const WorldVector<double>& x) const 
  {
    double normx = sqrt(x * x); // ||x||
    WorldMatrix<double> J;
    J.vecProduct(x, x); // xXx
    J *= - 1.0 / normx / normx; // - xXx/||x||^2
    for (int i = 0; i < 3; i++) J[i][i] += 1.0; // I - xXx/||x||^2
    return J / normx ;  // (I - xXx/||x||^2) / ||x||
  }


};




// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());
  
  DofEdgeVector alpha(edgeMesh, "alpha");
  alpha.set(new Alpha_d());

  //DofEdgeVector lb_alpha(edgeMesh, "lb_alpha");
  //alpha.set(new LbAlpha_d());
  //alpha.writeFile("output/lbalpha.vtu");


  DofEdgeVector dxyz(edgeMesh, "dxyz");
  dxyz.set(new DXYZ_d());

  //DofEdgeVector lcb_dxyz(edgeMesh, "lcb_dxyz");
  //lcb_dxyz.set(new LcbDXYZ_d());
  //lcb_dxyz.writeFile("output/lcbdxyz.vtu");

  DofEdgeVector initSol = dxyz + alpha;
  initSol.setName("initSol");
  initSol.writeFile("output/initSol.vtu");
  DOFVector< WorldVector<double> > initSolSharp = initSol.getSharpFaceAverage();
  io::VtkVectorWriter::writeFile(initSolSharp, "output/initSolSharp.vtu");


  DecProblemStat decSphere(&sphere, edgeMesh);

  DecProblemInstat sphereInstat(&decSphere);

  double minusOne = -1.0;
  // -Beltrami
  EdgeOperator Beltrami;
  Beltrami.addTerm(new LaplaceBeltramiAtEdges());
  decSphere.addMatrixOperator(Beltrami, 0, 0, &minusOne);
  // -CoBeltrami
  EdgeOperator CoBeltrami;
  CoBeltrami.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(CoBeltrami, 1, 1, &minusOne);
  // DeRham
  //EdgeOperator DeRham;
  //DeRham.addTerm(new LaplaceBeltramiAtEdges(-1.));
  //DeRham.addTerm(new LaplaceCoBeltramiAtEdges(-1.));
  //decSphere.addMatrixOperator(DeRham, 2, 2);
  decSphere.addMatrixOperator(Beltrami, 2, 2, &minusOne);
  decSphere.addMatrixOperator(CoBeltrami, 2, 2, &minusOne);

  // time derivatives approx
  EdgeOperator I0Operator;
  I0Operator.addTerm(new IdentityAtEdges());
  I0Operator.setUhOld(initSol);
  decSphere.addMatrixOperator(I0Operator, 0, 0,sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(I0Operator, 0,sphereInstat.getInvTauPtr());

  EdgeOperator I1Operator;
  I1Operator.addTerm(new IdentityAtEdges());
  I1Operator.setUhOld(initSol);
  decSphere.addMatrixOperator(I1Operator, 1, 1,sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(I1Operator, 1,sphereInstat.getInvTauPtr());

  EdgeOperator I2Operator;
  I2Operator.addTerm(new IdentityAtEdges());
  I2Operator.setUhOld(initSol);
  decSphere.addMatrixOperator(I2Operator, 2, 2,sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(I2Operator, 2,sphereInstat.getInvTauPtr());


  //decSphere.assembleSystem();
  //SparseMatrix sysMat = decSphere.getSysMat();
  //DenseVector rhs = decSphere.getRhs();
  //cout  << size(rhs) << endl;
  //int n = sysMat.num_rows();
  //double eps = 0.001;
  //for (int i = 0; i < n; ++i) {
  //  for (int j = 0; j < n; ++j) {
  //    if (i == j && sysMat[i][j] != 10) cout <<  sysMat[i][j] << endl;
  //    if (i != j && sysMat[i][j] != 0)  cout <<  sysMat[i][j] << endl;
  //  }
  //  if (rhs[i] != 10) cout << i << " : " <<  rhs[i] << endl;
  //}
  ////cout << sysMat << endl;

  //decSphere.solveDeprecated();
  



  sphereInstat.solve();
  

  //decSphere.writeSolution();


  //SparseMatrix sysMat = decSphere.getSysMat();
  //int n = edgeMesh->getNumberOfEdges();
  //mtl::dense2D<double> dmat(3*n, 3*n);
  //dmat = sysMat;
  //DenseVector eig(3*n);
  //eig = mtl::matrix::qr_algo(dmat,10);
  //cout  << "EVs: " << eig << endl;
  //mtl::dense2D<double> dmat(3*n, 3*n);
  //dmat = sysMat;
  //for (int i = 0; i < 3; ++i) {
  //  for (int j = 0; j < 3; ++j) {
  //    cout << i << "," << j << " Block ***********************" << endl;
  //    cout << sub_matrix(dmat, i*n, (i+1)*n, j*n, (j+1)*n) << endl;
  //  }
  //}
  //cout << sysMat << endl;

  //DenseVector rhs = decSphere.getRhs();
  //for (int i = 0; i < 3; ++i) {
  //  cout << i << "th subvector *****************" << endl;
  //  cout << sub_vector(rhs, i*n, (i+1)*n) << endl;
  //}
  //cout << rhs << endl;

  //DenseVector fsol = decSphere.getFullSolution();
  //cout << fsol << endl;



  //// === create adapt info ===
  //AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  //// === create adapt ===
  //AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
	//				       &sphere,
	//				       adaptInfo);
  

  // ===== start adaption loop =====
  //adapt->adapt();

  //sphere.writeFiles(adaptInfo, true);
  
  AMDiS::finalize();
}


