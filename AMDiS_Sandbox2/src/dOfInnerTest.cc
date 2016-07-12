#include "Dec.h"
#include "SphereProjection.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// <df,[p,q]>
class DZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = z
    return  q[2] - p[2];
  }
};


class DNorm2DZ_Sphere : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > {
  public:
  DNorm2DZ_Sphere() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const {
    // d(||dz||^2) = d(1 - z^2) = -d(z^2)
    return p[2]*p[2] - q[2]*q[2];
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

// Definition of df //
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dff = new DZ();
  
  DofEdgeVector df(edgeMesh, "df");
  df.set(dff);
  df.writeSharpFile("output/df.vtu", &sphere);

// Definition of d||df||^2 // 
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dn2dff = new DNorm2DZ_Sphere();
  
  DofEdgeVector dn2df(edgeMesh, "dOfNormSquareOfdf");
  dn2df.set(dn2dff);
  dn2df.writeSharpFile("output/dn2df.vtu", &sphere);
  dn2df.writeFile("output/dn2dfForm.vtu");




  DofEdgeVector dn2dfh = df.exteriorDerivativOfNorm2();
  dn2dfh.writeSharpFile("output/dn2dfh.vtu", &sphere);
  dn2dfh.writeFile("output/dn2dfhForm.vtu");

  //decSphere.assembleSystem();

  //decSphere.solve();

  //decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  double errL2Rel   = dn2df.errorL2Rel (dn2dfh);
  double errMaxRel  = dn2df.errorMaxRel(dn2dfh);

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  cout << endl;
  cout << "rotf: RelError L2:  " << errL2Rel << endl;
  cout << "rotf: RelError Max: " << errMaxRel << endl;


  AMDiS::finalize();
}
