#include "Dec.h"
#include "SphereProjection.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// <df,[p,q]>
class Df : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Df() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = z
    return  q[2] - p[2];
    //return  q[0]*q[2] - p[0]*p[2];
  }
};

class Rotz : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Rotz() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> rotf;
    rotf[0] = y;
    rotf[1] = -x;
    rotf[2] = 0.0;
    return  rotf * vec;
  }
};

class RotzEllipt : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotzEllipt() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> rotf;
    rotf[0] = (-12.*y*sqrt((4.*std::pow(x,4)*std::pow(z,2) + 4.*std::pow(y,2)*(9. + (-18. + std::pow(y,2))*std::pow(z,2) + 9.*std::pow(z,4)) + std::pow(x,2)*(9. + 2.*(-9. + 4.*std::pow(y,2))*std::pow(z,2) + 9.*std::pow(z,4)))/std::pow(-1. + std::pow(z,2),2)))/(-45. + 27.*std::pow(x,2) - 27.*std::pow(y,2) + 37.*std::pow(z,2));
    rotf[1] = (6.*x*sqrt((4.*std::pow(x,4)*std::pow(z,2) + 4.*std::pow(y,2)*(9. + (-18. + std::pow(y,2))*std::pow(z,2) + 9.*std::pow(z,4)) + std::pow(x,2)*(9. + 2.*(-9. + 4.*std::pow(y,2))*std::pow(z,2) + 9.*std::pow(z,4)))/std::pow(-1. + std::pow(z,2),2)))/(-45. + 27.*std::pow(x,2) - 27.*std::pow(y,2) + 37.*std::pow(z,2));
    rotf[2] = 0.0;
    return  rotf * vec;
  }
};

// x * z
class f : public AbstractFunction<double, WorldVector<double> > {
  public:
  f() : AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& v) const {
    return v[2];
    //return v[0]*v[2];
  }
};


int main(int argc, char* argv[])
{
  FUNCNAME("main");
  
  AMDiS::init(argc, argv);
  
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  SphereProject sproj(42, VOLUME_PROJECTION);
  
  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  //DecProblemStat decSphere(&sphere, edgeMesh);

// Definition of rotf //
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dff = new Df();
  
  DofEdgeVector rotf(edgeMesh, "rotf");
  rotf.setDual(dff);
  //rotf.interpolGL4(new Rotz(), sproj.getProjection(), sproj.getJProjection());
  //rotf.interpolGL4(new RotzEllipt(), sproj.getProjection(), sproj.getJProjection());
  rotf.writeSharpFile("output/rot.vtu", &sphere);

// Definition of f // 
  
  AbstractFunction<double, WorldVector<double> > *ff = new f();
  DofVertexVector fvec(edgeMesh, "f");
  fvec.interpol(ff);

  DofEdgeVector rotfh = *(fvec.rotOnEdges());
  rotfh.writeSharpFile("output/roth.vtu", &sphere);

  //decSphere.assembleSystem();

  //decSphere.solve();

  //decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  double errL2Rel_p = rotf.errorL2Rel(rotfh);
  double errMaxRel_p = rotf.errorMaxRel(rotfh);

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  cout << endl;
  cout << "rotf: RelError L2:  " << errL2Rel_p << endl;
  cout << "rotf: RelError Max: " << errMaxRel_p << endl;


  AMDiS::finalize();
}
