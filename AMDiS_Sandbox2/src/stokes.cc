#include "Dec.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// <df,[p,q]>
class DXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = z
    //return  q[2] - p[2];
    return q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
  }
};

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

class GaussCurv_Sphere : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_Sphere() : AbstractFunction<double, EdgeElement >(){}

  double operator()(const EdgeElement& eel) const {
    return 1.0;
  }
};

class GaussCurv_Ellipsoid : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_Ellipsoid() : AbstractFunction<double, EdgeElement >(){}

  double operator()(const EdgeElement& eel) const {
    WorldVector<double> x = eel.infoLeft->getEdgeCenter(eel.dofEdge);
    double nSqrt= 81. + 972.*x[1]*x[1] -20.*x[2]*x[2];
    return 11664./nSqrt/nSqrt;
  }
};


int main(int argc, char* argv[])
{
  FUNCNAME("main");
  
  AMDiS::init(argc, argv);
  
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);
  
  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);

// gaussCurv on edge circumcenter

  //AbstractFunction<double, EdgeElement > *Kf = new GaussCurv_Sphere();
  AbstractFunction<double, EdgeElement > *Kf = new GaussCurv_Ellipsoid();
  DofEdgeVector K(edgeMesh, "K");
  K.set(Kf);

// Definition of rhs //
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dzf = new DZ();
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dxyzf = new DXYZ();
  
  DofEdgeVector FPot(edgeMesh, "FPot");
  FPot.set(dzf);
  FPot.writeSharpFile("output/FPot.vtu", &sphere);

  DofEdgeVector FRot(edgeMesh, "FRot");
  FRot.setDual(dxyzf);
  FRot.writeSharpFile("output/FRot.vtu", &sphere);

  double mu = 1.0;

  EdgeOperator muRotrotu;
  muRotrotu.addTerm(new LaplaceBeltramiAtEdges(mu));
  decSphere.addMatrixOperator(muRotrotu,0,0);

  EdgeOperator muKu;
  //muKu.addTerm(new IdentityAtEdges(mu));
  muKu.addTerm(new EdgeVecAtEdges(&K, NULL, 2.0*mu));
  decSphere.addMatrixOperator(muKu,0,0);

  EdgeOperator dp;
  dp.addTerm(new ExteriorDerivativeAtEdges(-1.0));
  decSphere.addMatrixOperator(dp,0,1);

  EdgeOperator rhs0;
  rhs0.addTerm(new EdgeVecAtEdges(&FPot));
  rhs0.addTerm(new EdgeVecAtEdges(&FRot));
  decSphere.addVectorOperator(rhs0,0);

  VertexOperator divu;
  divu.addTerm(new DivAtVertices());
  decSphere.addMatrixOperator(divu,1,0);



  decSphere.assembleSystem();
  //cout << decSphere.getSysMat() << endl;
  //cout << decSphere.getRhs() << endl;


  decSphere.setValAtDof(1, 0, 0.);
  //cout << decSphere.getSysMat() << endl;
  //cout << decSphere.getRhs() << endl;

  decSphere.solve();

  decSphere.writeSolution();
  
  


  AMDiS::finalize();
}
