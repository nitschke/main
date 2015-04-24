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
  //alpha.writeFile("output/alpha.vtu");

  //DofEdgeVector lb_alpha(edgeMesh, "lb_alpha");
  //alpha.set(new LbAlpha_d());
  //alpha.writeFile("output/lbalpha.vtu");


  DofEdgeVector dxyz(edgeMesh, "dxyz");
  dxyz.set(new DXYZ_d());
  //dxyz.writeFile("output/dxyz.vtu");

  //DofEdgeVector lcb_dxyz(edgeMesh, "lcb_dxyz");
  //lcb_dxyz.set(new LcbDXYZ_d());
  //lcb_dxyz.writeFile("output/lcbdxyz.vtu");

  DofEdgeVector initSol = dxyz + alpha;

  DecProblemStat decSphere(&sphere, edgeMesh);

  DecProblemInstat sphereInstat(&decSphere);

 // LaplaceDeRham
  double minusOne = -1.0;
  EdgeOperator LaplaceOperator;
  LaplaceOperator.addTerm(new LaplaceBeltramiAtEdges());
  //LaplaceOperator.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceOperator, 0, 0, &minusOne);

  // time derivative approx
  EdgeOperator IOperator;
  IOperator.addTerm(new IdentityAtEdges());
  IOperator.setUhOld(&initSol);
  decSphere.addMatrixOperator(IOperator, 0, 0,sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(IOperator, 0,sphereInstat.getInvTauPtr());




  sphereInstat.solve();
  
  decSphere.writeSolution();


  //SparseMatrix sysMat = decSphere.getSysMat();
  //cout << sysMat << endl;

  //DenseVector rhs = decSphere.getRhs();
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


