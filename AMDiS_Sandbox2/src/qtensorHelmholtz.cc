#include "Dec.h"
#include "SphereProjection.h"


using namespace std;
using namespace AMDiS;
using namespace dec;

class DZ2 : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZ2() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[2]*q[2] - p[2]*p[2];
  }
};

class N2q : public AbstractFunction<double, WorldVector<double> > {
public:
  N2q() :  AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& coords) const {
    double t = 1.0 - coords[2]*coords[2];
    return 8.0 * t * t;
  }
};

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  SphereProject proj(42, VOLUME_PROJECTION);


  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);

  
  //Laplace
  EdgeOperator LaplaceB;
  LaplaceB.addTerm(new LaplaceBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceB, 0, 0);
  decSphere.addMatrixOperator(LaplaceB, 1, 1);

  EdgeOperator LaplaceCB;
  LaplaceCB.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceCB, 0, 0);
  decSphere.addMatrixOperator(LaplaceCB, 1, 1);

  //Ids
  EdgeOperator Id00;
  Id00.addTerm(new IdentityAtEdges(2.0));
  decSphere.addMatrixOperator(Id00,0,0);

  EdgeOperator Id11;
  Id11.addTerm(new IdentityAtEdges(7.0));
  decSphere.addMatrixOperator(Id11,1,1);

  EdgeOperator Id01;
  Id01.addTerm(new IdentityAtEdges(-1.0));
  decSphere.addMatrixOperator(Id01,0,1);

  //RHS
  DofEdgeVector dz2(edgeMesh, "dz2");
  dz2.set(new DZ2());
  dz2.writeSharpFile("output/dz2.vtu", &sphere);
  EdgeOperator RHS1;
  RHS1.addTerm(new EdgeVecAtEdges(&dz2,-4.0));
  decSphere.addVectorOperator(RHS1,1);

  //decSphere.setValAtDof(0, 0, 500, 500, 0.0);
  //decSphere.setValAtDof(0, 0, 250, 250, 0.0);
  //decSphere.setValAtDof(0, 0, 0, 0, 0.0);

  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();

  DofEdgeVector alpha = decSphere.getEdgeSolution(0);
  DofEdgeVector halpha = alpha.hodgeDual();
  DofEdgeVector beta = decSphere.getEdgeSolution(1);
  DofEdgeVector hbeta = beta.hodgeDual();
  DofEdgeVector hdz2 = dz2.hodgeDual();

  double errL2 = (256./15.)*M_PI - 2.*DofEdgeVectorPD::L2Inner(alpha, halpha, beta, hbeta) - 16.*DofEdgeVectorPD::L2Inner(alpha, halpha, dz2, hdz2);
  cout  << "L2-Error of q=D(alpha): " << errL2 << endl;

  DofVertexVector n2alpha = DofEdgeVectorPD::normOnVertices(alpha, halpha);
  DofVertexVector divalpha = alpha.divergence();
  DofVertexVector rotalpha = -1.0 * halpha.divergence();
  DofVertexVector n2qh = 2.0 * (*(n2alpha.laplace())
                                + 2.0 * n2alpha 
                                - 2.0 * DofEdgeVectorPD::interiorProdOnVertices(alpha, halpha, beta, hbeta)
                                - (divalpha * divalpha)
                                - (rotalpha * rotalpha));
  DofVertexVector n2q(edgeMesh, "n2q");
  n2q.interpol(new N2q());
  double norm2ErrMax = n2qh.errorMax(n2q);
  double norm2ErrL2 = n2qh.errorL2(n2q);
  cout << "norm2ErrMax: " << norm2ErrMax << endl;
  cout << "norm2ErrL2 : " << norm2ErrL2 << endl;
  
  DofVertexVector diffn2q = n2q - n2qh;
  diffn2q.writeFile("output/diffn2q.vtu");


  AMDiS::finalize();
}
