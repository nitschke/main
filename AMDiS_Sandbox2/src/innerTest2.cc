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

class G : public AbstractFunction<double, WorldVector<double> > {
  public:
  G() : AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& coords) const {
    return 1.0 - coords[2]*coords[2] + coords[2];
  }
};

class f : public AbstractFunction<double, WorldVector<double> > {
  public:
  f() : AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& v) const {
    return v[2];
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

  DecProblemStat decSphere(&sphere, edgeMesh);

// Definition of df //
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dff = new DZ();
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *rotff = new RotZ_Sphere();
  
  DofEdgeVector dfP(edgeMesh, "df");
  dfP.set(dff);
  DofEdgeVector dfD(edgeMesh, "Rotf");
  dfD.interpolGL4(rotff, sproj.getProjection(), sproj.getJProjection());
  DofEdgeVectorPD df(dfP, dfD);
  df.writeSharpFile("output/df.vtu", &sphere);


  // gamma = df
  // f = z

  // <df,gamma>
  VertexOperator interP;
  interP.addTerm(new InterProdPartAtVertices(&dfP));
  decSphere.addMatrixOperator(interP,0,1);

  VertexOperator interD;
  interD.addTerm(new InterProdPartAtVertices(&dfD));
  decSphere.addMatrixOperator(interD,0,2);

  // f
  VertexOperator fOp;
  fOp.addTerm(new IdentityAtVertices());
  decSphere.addMatrixOperator(fOp,0,0);

  // g = <df,gamma> + f = 1- z^2 +z (on Sphere)
  VertexOperator gOp;
  AbstractFunction<double, WorldVector<double> > *gf = new G();
  DofVertexVector gvec(edgeMesh, "f");
  gvec.interpol(gf);
  gOp.addTerm(new VertexVecAtVertices(&gvec));
  decSphere.addVectorOperator(gOp,0);

  // df = df
  EdgeOperator dfPOp;
  dfPOp.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(dfPOp,1,1);
  
  EdgeOperator rhsPOp;
  rhsPOp.addTerm(new EdgeVecAtEdges(&dfP));
  decSphere.addVectorOperator(rhsPOp,1);

  // *df = *df
  EdgeOperator dfDOp;
  dfDOp.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(dfDOp,2,2);
  
  EdgeOperator rhsDOp;
  rhsDOp.addTerm(new EdgeVecAtEdges(&dfD));
  decSphere.addVectorOperator(rhsDOp,2);

  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  double errL2Rel_P = decSphere.getEdgeSolution(1).errorL2Rel(dfP);
  double errMaxRel_P = decSphere.getEdgeSolution(1).errorMaxRel(dfP);

  double errL2Rel_D = decSphere.getEdgeSolution(2).errorL2Rel(dfD);
  double errMaxRel_D = decSphere.getEdgeSolution(2).errorMaxRel(dfD);

  cout << endl;
  cout << "P: RelError L2:  " << errL2Rel_P << endl;
  cout << "P: RelError Max: " << errMaxRel_P << endl;
  cout << endl;
  cout << "D: RelError L2:  " << errL2Rel_D << endl;
  cout << "D: RelError Max: " << errMaxRel_D << endl;

  DofVertexVector fhvec = decSphere.getVertexSolution(0);
  AbstractFunction<double, WorldVector<double> > *ff = new f();
  DofVertexVector fvec(edgeMesh, "f");
  fvec.interpol(ff);
  double errMaxRel_f = fhvec.errorMaxRel(fvec);
  double errL2Rel_f = fhvec.errorL2Rel(fvec);
  cout << endl;
  cout << "f: RelError L2:  " << errL2Rel_f << endl;
  cout << "f: RelError Max: " << errMaxRel_f << endl;  


  AMDiS::finalize();
}
