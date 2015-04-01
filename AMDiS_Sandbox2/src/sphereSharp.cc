#include "AMDiS.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"
#include "io/VtkVectorWriter.h"
#include "io/ElementFileWriter.h"

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
  //cout << *edgeMesh << endl;

 
  cout << endl << "********** Gamma *************" << endl;
  DofEdgeVector gammad(edgeMesh, "gamma_d");
  gammad.interpolGL4(new Gamma(), new Proj(), new JProj());
  DofEdgeVector gammadExact(edgeMesh, "gamma_d_exact");
  gammadExact.set(new Gamma_d());
  DofEdgeVector diff = gammad - gammadExact;
  cout << "Error GL4: " << diff.l2Norm() << endl;
  gammad.interpolLinTrapz(new Gamma());
  diff = gammad - gammadExact;
  cout << "Error LinTrapz: " << diff.l2Norm() << endl;
  gammad.interpolLinMidpoint(new Gamma());
  diff = gammad - gammadExact;
  cout << "Error LinMidpoint: " << diff.l2Norm() << endl;
  gammad.interpolMidpoint(new Gamma(), new Proj());
  diff = gammad - gammadExact;
  cout << "Error Midpoint: " << diff.l2Norm() << endl;
  gammad.interpolNC(new Gamma(), 3);
  diff = gammad - gammadExact;
  cout << "Error LinSimpson: " << diff.l2Norm() << endl;
  gammad.interpolNC(new Gamma(), 3, new Proj());
  diff = gammad - gammadExact;
  cout << "Error Simpson: " << diff.l2Norm() << endl;
  gammad.interpolNC(new Gamma(), 7);
  diff = gammad - gammadExact;
  cout << "Error LinWeddle: " << diff.l2Norm() << endl;
  gammad.interpolNC(new Gamma(), 7, new Proj());
  diff = gammad - gammadExact;
  cout << "Error Weddle: " << diff.l2Norm() << endl;

  //DOFVector< WorldVector<double> > gammaSharp = gammadExact.getSharpEdgeRingLinMod();
  DOFVector< WorldVector<double> > gammaSharp = gammadExact.getSharpFaceAverage();
  AMDiS::io::VtkVectorWriter::writeFile(gammaSharp, string("output/gammaSharp.vtu"));

  cout << endl << "********** Alpha *************" << endl;
  DofEdgeVector alphadGL4(edgeMesh, "alpha_d_GL4");
  alphadGL4.interpolGL4(new Alpha(), new Proj(), new JProj());
  
  DofEdgeVector alphaExact(edgeMesh, "alphaExact");
  alphaExact.set(new Alpha_d());
  DofEdgeVector diff2 = alphadGL4 - alphaExact;
  cout << "Error GL4: " << diff2.l2Norm() << endl;

  DofEdgeVector alphad(edgeMesh, "alpha_d");
  alphad.interpolLinTrapz(new Alpha());
  diff = alphad - alphadGL4;
  cout << "Error LinTrapz: " << diff.l2Norm() << endl;
  alphad.interpolLinMidpoint(new Alpha());
  diff = alphad - alphadGL4;
  cout << "Error LinMidpoint: " << diff.l2Norm() << endl;
  alphad.interpolMidpoint(new Alpha(), new Proj());
  diff = alphad - alphadGL4;
  cout << "Error Midpoint: " << diff.l2Norm() << endl;
  alphad.interpolNC(new Alpha(), 3);
  diff = alphad - alphadGL4;
  cout << "Error LinSimpson: " << diff.l2Norm() << endl;
  alphad.interpolNC(new Alpha(), 3, new Proj());
  diff = alphad - alphadGL4;
  cout << "Error Simpson: " << diff.l2Norm() << endl;
  alphad.interpolNC(new Alpha(), 7);
  diff = alphad - alphadGL4;
  cout << "Error LinWeddle: " << diff.l2Norm() << endl;
  alphad.interpolNC(new Alpha(), 7, new Proj());
  diff = alphad - alphadGL4;
  cout << "Error Weddle: " << diff.l2Norm() << endl;

  //DOFVector< WorldVector<double> > alphaSharp = alphadGL4.getSharpFaceAverage();
  //DOFVector< WorldVector<double> > alphaSharp = alphadGL4.getSharpEdgeRingLinMod();
  //DOFVector< WorldVector<double> > alphaSharp = alphadGL4.getSharpHirani();
  DOFVector< WorldVector<double> > alphaSharp = alphadGL4.getSharpEdgeRing();
  AMDiS::io::VtkVectorWriter::writeFile(alphaSharp, string("output/alphaSharp.vtu"));

  //map<int, std::vector<double> > alphaFaceSharp = alphadGL4.getSharpOnFaces();
  //AMDiS::io::ElementFileWriter::writeFile(alphaFaceSharp, sphere.getFeSpace()->getMesh(), "output/alphaFaceSharp");

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


