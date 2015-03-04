#include "AMDiS.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"
#include "io/VtkVectorWriter.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

// 1-form
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

// exact 1-form : beta = dz
class Beta : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Beta() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    WorldVector<double> gradF;
    gradF[0] = - x[0]*x[2];
    gradF[1] = - x[1]*x[2];
    gradF[2] = 1.0 - x[2]*x[2];
    return  gradF * vec;
  }
};

// <beta,[p,q]>
class Beta_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Beta_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[2] - p[2];
  }
};

// exact 1-form : gamma = d(xyz)
class Gamma : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Gamma() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

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

// <gamma,[p,q]>
class Gamma_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Gamma_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
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

  const EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  //DofEdgeVector alphad(edgeMesh, "alpha_d");
  //alphad.interpolGL4(new Alpha(), new Proj(), new JProj());
  
  //DofEdgeVector betad(edgeMesh, "beta_d");
  //betad.interpolGL4(new Beta(), new Proj(), new JProj());
  //DofEdgeVector betadExact(edgeMesh, "beta_d_exact");
  //betadExact.set(new Beta_d());
  //DofEdgeVector diff = betad - betadExact;
  //cout << "Error: " << diff.L2Norm() << endl;
  
  DofEdgeVector gammad(edgeMesh, "gamma_d");
  gammad.interpolGL4(new Gamma(), new Proj(), new JProj());
  DofEdgeVector gammadExact(edgeMesh, "gamma_d_exact");
  gammadExact.set(new Gamma_d());
  DofEdgeVector diff = gammad - gammadExact;
  cout << "Error GL4: " << diff.L2Norm() << endl;
  gammad.interpolSimple(new Gamma());
  diff = gammad - gammadExact;
  cout << "Error Simple: " << diff.L2Norm() << endl;



  //// === create adapt info ===
  //AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  //// === create adapt ===
  //AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
	//				       &sphere,
	//				       adaptInfo);
  

  // ===== start adaption loop =====
  //adapt->adapt();

  //sphere.writeFiles(adaptInfo, true);
  


  //DOFVector<WorldVector<double> > WDV(sphere.getFeSpace(),"vector");
  //WDV.interpol(new Alpha());
  //AMDiS::io::VtkVectorWriter::writeFile(WDV, string("output/vector.vtu"));

  AMDiS::finalize();
}


