#include "Dec.h"

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
  }
};

// <dLaplace(f),[p,q]>
class DLBfSphere : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DLBfSphere() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // Laplace(f)(X) = -2z
    return  -2.*(q[2] - p[2]);
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

// Definition of p //
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dff = new Df();
  DofEdgeVector df(edgeMesh, "df");
  df.set(dff);

  DofEdgeVector stardf(edgeMesh, "stardf");
  stardf.setDual(dff);
  
  //p = df + *df
  DofEdgeVector p(edgeMesh, "p");
  p.set(0.0);
  p = df + stardf;

// Definition of Laplace(p) - p  = -LaplaceDeRham(p) -p

  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dLBff = new DLBfSphere();
  DofEdgeVector dLBf(edgeMesh, "dLBf");
  dLBf.set(dLBff);

  DofEdgeVector stardLBf(edgeMesh, "stardLBf");
  stardLBf.setDual(dLBff);
  
  //laplace(p) - p = dLaplace(f) + *dLaplace(f) - p =  dLaplace(f) + *dLaplace(f) - ( df + *df )
  DofEdgeVector rhs(edgeMesh, "rhs");
  rhs.set(0.0);
  rhs = dLBf + stardLBf - p;


 // -LaplaceDeRham
  EdgeOperator LaplaceOperator;
  LaplaceOperator.addTerm(new LaplaceBeltramiAtEdges());
  LaplaceOperator.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceOperator, 0, 0);

  EdgeOperator IOperator;
  IOperator.addTerm(new IdentityAtEdges(-1.0));
  decSphere.addMatrixOperator(IOperator, 0, 0);


  EdgeOperator VOperator;
  VOperator.addTerm(new EdgeVecAtEdges(&rhs));
  decSphere.addVectorOperator(VOperator, 0);

  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  double errL2Rel = decSphere.getSolution().errorL2Rel(p);
  double errMaxRel = decSphere.getSolution().errorMaxRel(p);

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  cout << endl;
  cout << "RelError L2:  " << errL2Rel << endl;
  cout << "RelError Max: " << errMaxRel << endl;

  cout << endl;
  cout << "maxDiameter,maxLength,errorL2Rel,errorMaxRel" << endl;
  cout << hdia << "," << hlen << "," << errL2Rel << ", " << errMaxRel << endl;

  AMDiS::finalize();
}
