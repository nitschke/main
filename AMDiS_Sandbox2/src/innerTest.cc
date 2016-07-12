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

class RotZ_Sphere : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotZ_Sphere() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    double con2 = 1.0; //contra coords
    WorldVector<double> conBasis2; //contra basis vectors
    conBasis2[0] = x[1];
    conBasis2[1] = -x[0];
    conBasis2[2] = 0.0;
    return (con2 * conBasis2) * vec;
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

class Norm2DZ : public AbstractFunction<double, WorldVector<double> > {
  public:
  Norm2DZ() : AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& coords) const {
    return 1.0 - coords[2]*coords[2];
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
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *rotff = new RotZ_Sphere();
  
  DofEdgeVector dfP(edgeMesh, "df");
  dfP.set(dff);
  DofEdgeVector dfD(edgeMesh, "Rotf");
  dfD.interpolGL4(rotff, sproj.getProjection(), sproj.getJProjection());
  DofEdgeVectorPD df(dfP, dfD);
  df.writeSharpFile("output/df.vtu", &sphere);

// Definition of ||df||^2

  AbstractFunction<double, WorldVector<double> > *norm2dff = new Norm2DZ();

  DofVertexVector norm2df(edgeMesh, "norm2df");
  norm2df.interpol(norm2dff);
  norm2df.writeFile("output/norm2df.vtu");


// calc norm

  DofVertexVector norm2dfh = df.interiorProdOnVertices(df);
  norm2dfh.writeFile("output/norm2dfh.vtu");


// Definition of d||df||^2 // 
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dn2dff = new DNorm2DZ_Sphere();
  
  DofEdgeVector dn2df(edgeMesh, "dOfNormSquareOfdf");
  dn2df.set(dn2dff);
  dn2df.writeSharpFile("output/dn2df.vtu", &sphere);
  dn2df.writeFile("output/dn2dfForm.vtu");

  DofEdgeVector dn2dfh(*(norm2dfh.exteriorDerivative()));
  dn2dfh.writeSharpFile("output/dn2dfh.vtu", &sphere);
  dn2dfh.writeFile("output/dn2dfhForm.vtu");

  //decSphere.assembleSystem();

  //decSphere.solve();

  //decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  double errL2Rel_norm2   = norm2dfh.errorL2Rel (norm2df);
  double errMaxRel_norm2  = norm2dfh.errorMaxRel(norm2df);
  double errL2Rel_dnorm2   = dn2df.errorL2Rel (dn2dfh);
  double errMaxRel_dnorm2  = dn2df.errorMaxRel(dn2dfh);

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  cout << endl;
  cout << "norm2: RelError L2:  " << errL2Rel_norm2 << endl;
  cout << "norm2: RelError Max: " << errMaxRel_norm2 << endl;

  cout << endl;
  cout << "d(norm2): RelError L2:  " << errL2Rel_dnorm2 << endl;
  cout << "d(norm2): RelError Max: " << errMaxRel_dnorm2 << endl;


  AMDiS::finalize();
}
