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
  
  DofEdgeVector p(edgeMesh, "p");
  p.set(dff);

  DofEdgeVector starp(edgeMesh, "stardf");
  starp.setDual(dff);

  DofEdgeVectorPD::writeSharpOnEdgeFile(p, starp, "output/pOnEdges.vtu");

  // rot(p) = -div(*p) // rot(df) = *ddf = 0
  DofEdgeVector rotp = -1.0 * starp.divOnEdgeCenter();
  DofEdgeVector rotpExact(edgeMesh,"rotpEx");
  rotpExact.set(0.0);
  cout  << endl;
  cout  << "************ rot = -div* **************************" << endl;
  cout << "RelError L2:  " << rotp.errorL2(rotpExact) << endl;
  cout << "RelError Max: " << rotp.errorMax(rotpExact) << endl;
  cout  << "***************************************************" << endl;
  

  //EdgeOperator VOperator;
  //VOperator.addTerm(new EdgeVecAtEdges(&rhs));
  //decSphere.addVectorOperator(VOperator, 0);

  //decSphere.assembleSystem();

  //decSphere.solve();

  //decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  //double errL2Rel = decSphere.getSolution().errorL2Rel(p);
  //double errMaxRel = decSphere.getSolution().errorMaxRel(p);

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  //cout << endl;
  //cout << "RelError L2:  " << errL2Rel << endl;
  //cout << "RelError Max: " << errMaxRel << endl;

  //cout << endl;
  //cout << "maxDiameter,maxLength,errorL2Rel,errorMaxRel" << endl;
  //cout << hdia << "," << hlen << "," << errL2Rel << ", " << errMaxRel << endl;

  AMDiS::finalize();
}
