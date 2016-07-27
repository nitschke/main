#include "Dec.h"
#include "phiProjection.h"
#include "EllipsoidProjection.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

class Df : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Df() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = z
    //return  q[2] - p[2];
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
    //return  q[0]*q[2] - p[0]*p[2];
  }
};

class GradfEllipt : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  GradfEllipt() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> dfvec;
    ////Sphere dz
    //dfvec[0] = -1.*x*z;
    //dfvec[1] = -1.*y*z;
    //dfvec[2] = 1. - 1.*std::pow(z,2);
    //// Ellipsoid dz
    //dfvec[0] = (72.*x*z)/(-405. + 243.*std::pow(x,2) - 972.*std::pow(y,2) + 148.*std::pow(z,2));
    //dfvec[1] = (288.*y*z)/(-405. + 243.*std::pow(x,2) - 972.*std::pow(y,2) + 148.*std::pow(z,2));
    //dfvec[2] = (9.*(-45. + 27.*std::pow(x,2) - 108.*std::pow(y,2) + 20.*std::pow(z,2)))/(-405. + 243.*std::pow(x,2) - 972.*std::pow(y,2) + 148.*std::pow(z,2));
    ////sphere d(xyz)
    //for (int i = 0; i < 3; i++) {
    //  int ii = (i+1)%3;
    //  int iii = (i+2)%3;
    //  dfvec[i] = coords[ii] * coords[iii] * (1.0 - 3.0*coords[i]*coords[i]);
    //}
    // Ellipsoid d(xyz)
    dfvec[0] = (2.*y*z*(153. - 513.*std::pow(x,2) + 684.*std::pow(y,2) - 52.*std::pow(z,2)))/(405. - 243.*std::pow(x,2) + 972.*std::pow(y,2) - 148.*std::pow(z,2));
    dfvec[1] = (-1.*x*z*(63. + 99.*std::pow(x,2) - 1188.*std::pow(y,2) + 4.*std::pow(z,2)))/(-405. + 243.*std::pow(x,2) - 972.*std::pow(y,2) + 148.*std::pow(z,2));
    dfvec[2] = (27.*x*y*(-15. + 9.*std::pow(x,2) - 36.*std::pow(y,2) + 20.*std::pow(z,2)))/(-405. + 243.*std::pow(x,2) - 972.*std::pow(y,2) + 148.*std::pow(z,2));
    return  dfvec * vec;
  }
};

//// **********  Sphere ********************////////
// for phi-projection
class Phi : public AbstractFunction<double, WorldVector<double> >
{
public:
  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& x) const 
  {
    return 0.5 * (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1.0);
  }
};

class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}

  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> rval(x);
    return rval;
  }
};
////**********************************************/////////////

////// **********  Ellipsoid ********************////////
//// for phi-projection
//class Phi : public AbstractFunction<double, WorldVector<double> >
//{
//public:
//  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}
//
//  double operator()(const WorldVector<double>& x) const 
//  {
//    return 0.5 * (x[0]*x[0] + 4.0*x[1]*x[1] + (4.0/9.0)*x[2]*x[2] - 1.0);
//  }
//};
//
//class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
//{
//public:
//  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}
//
//  WorldVector<double> operator()(const WorldVector<double>& x) const 
//  {
//    WorldVector<double> rval(x);
//    rval[1] *= 4.0;
//    rval[2] *= 4.0/9.0;
//    return rval;
//  }
//};
//////**********************************************/////////////

int main(int argc, char* argv[])
{
  FUNCNAME("main");
  
  AMDiS::init(argc, argv);
  
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  //PhiProject proj(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-15);
  EllipsoidProject proj(1, VOLUME_PROJECTION, 1.0, 0.5, 1.5);
  
  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DofEdgeVector df(edgeMesh, "df");
  df.set(new Df());
  df.writeSharpFile("output/df.vtu",&sphere);

  DofEdgeVector dfGL4(edgeMesh, "dfGL4");
  //dfGL4.interpolGL4(new GradfEllipt(), proj.getProjection(), proj.getJProjection(1.E-4));
  dfGL4.interpolGL4(new GradfEllipt(), proj.getProjection(), proj.getJProjection());
  dfGL4.writeSharpFile("output/dfGL4.vtu",&sphere);

  DofEdgeVector dfLM(edgeMesh, "dfLM");
  dfLM.interpolLinMidpoint(new GradfEllipt());
  dfLM.writeSharpFile("output/dfLM.vtu",&sphere);

  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  double errL2Rel_GL4 = dfGL4.errorL2Rel(df);
  double errMaxRel_GL4 = dfGL4.errorMaxRel(df);
  double errL2Rel_LM = dfLM.errorL2Rel(df);
  double errMaxRel_LM = dfLM.errorMaxRel(df);

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  cout << endl;
  cout << "dfGL4: RelError L2:  " << errL2Rel_GL4 << endl;
  cout << "dfGL4: RelError Max: " << errMaxRel_GL4 << endl;

  cout << endl;
  cout << "dfLM: RelError L2:  " << errL2Rel_LM << endl;
  cout << "dfLM: RelError Max: " << errMaxRel_LM << endl;

  AMDiS::finalize();
}
