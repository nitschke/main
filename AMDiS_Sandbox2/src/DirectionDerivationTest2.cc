#include "Dec.h"
#include "SphereProjection.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// q = 0.5 dz^2 = z*dz
class Q : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Q() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  0.5 * (q[2]*q[2] - p[2]*p[2]);
  }
};

//*q = z*Rot(z)
class HQ_Sphere : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  HQ_Sphere() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    double con2 = 1.0; //contra coords
    WorldVector<double> conBasis2; //contra basis vectors
    conBasis2[0] = x[1];
    conBasis2[1] = -x[0];
    conBasis2[2] = 0.0;
    return x[2] * (con2 * conBasis2) * vec;
  }
};

// Rot(z)
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

class LBZ : public AbstractFunction<double, WorldVector<double> > {
  public:
  LBZ() : AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& coords) const {
    return -2.0*coords[2];
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

  DecProblemStat decSphere(&sphere, edgeMesh);

// Definition of alpha = [*dz, -dz]//
  DofEdgeVector alphaP(edgeMesh, "alphaPrimal");
  DofEdgeVector alphaD(edgeMesh, "alphaDual");
  alphaP.interpolGL4(new RotZ_Sphere(), sproj.getProjection(), sproj.getJProjection());
  alphaD.set(new DZ());
  alphaD *= -1.0;

// Definition of q = 0.5 [d(z^2) ,  *d(z^2) ]
  DofEdgeVector qP(edgeMesh, "qPrimal");
  DofEdgeVector qD(edgeMesh, "qDual");
  qP.set(new Q());
  qD.interpolGL4(new HQ_Sphere(), sproj.getProjection(), sproj.getJProjection());
  
// Components [alpha, *alpha , r = rot(alpha) = - div(*alpha) , f = ||alpha||^2 = <alpha, alphaTilde>]

// rot(alpha) (*alphaTilde) + 0.5 df = q on Edges

  // + alpha - alphaTilde 
  EdgeOperator alphaId;
  alphaId.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(alphaId, 0, 0);
  EdgeOperator alphaRHS;
  alphaRHS.addTerm(new EdgeVecAtEdges(&alphaP));
  decSphere.addVectorOperator(alphaRHS, 0);

  EdgeOperator hAlphaTildeRotAlpha;
  hAlphaTildeRotAlpha.addTerm(new RotAtEdgeCenterAndEdgeVecAtEdges(&alphaD));
  decSphere.addMatrixOperator(hAlphaTildeRotAlpha, 0, 0);
  //rHAlphaTilde.addTerm(new EdgeVecAtEdges(&qP,-2.0));
  //decSphere.addVectorOperator(rHAlphaTilde, 0);
  //EdgeOperator hAlphaTildeDivHAlpha;
  //hAlphaTildeDivHAlpha.addTerm(new DivAtEdgeCenterAndEdgeVecAtEdges(&alphaD, -1.0));
  //decSphere.addMatrixOperator(hAlphaTildeDivHAlpha, 0, 1);

  EdgeOperator halfDF;
  halfDF.addTerm(new ExteriorDerivativeAtEdges(0.5));
  decSphere.addMatrixOperator(halfDF, 0, 2);

  EdgeOperator qRHS;
  qRHS.addTerm(new EdgeVecAtEdges(&qP));
  decSphere.addVectorOperator(qRHS, 0);

// -Rot(alpha) alphaTilde + 0.5 Rot(f) = *q  (Hodge Dual Equation)
  // + *alpha - *alpha
  EdgeOperator halphaId;
  halphaId.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(halphaId, 1, 1);
  EdgeOperator halphaRHS;
  halphaRHS.addTerm(new EdgeVecAtEdges(&alphaD));
  decSphere.addVectorOperator(halphaRHS, 1);

  EdgeOperator alphaTildeRotAlpha;
  alphaTildeRotAlpha.addTerm(new RotAtEdgeCenterAndEdgeVecAtEdges(&alphaP, -1.0));
  decSphere.addMatrixOperator(alphaTildeRotAlpha, 1, 0);
  //rAlphaTilde.addTerm(new EdgeVecAtEdges(&qD,-2.0));
  //decSphere.addVectorOperator(rAlphaTilde, 1);
  //EdgeOperator alphaTildeDivHAlpha;
  //alphaTildeDivHAlpha.addTerm(new DivAtEdgeCenterAndEdgeVecAtEdges(&alphaP));
  //decSphere.addMatrixOperator(alphaTildeDivHAlpha, 1, 1);

  EdgeOperator halfRotF;
  halfRotF.addTerm(new RotAtEdges(0.5));
  decSphere.addMatrixOperator(halfRotF, 1, 2);

  EdgeOperator hqRHS;
  hqRHS.addTerm(new EdgeVecAtEdges(&qD));
  decSphere.addVectorOperator(hqRHS, 1);


// j_alphaTildeSharp(alpha) + j_*alphaTildeSharp(*alpha) - f = 0    ( ---> <alpha,alphaTilde> - f = 0)
  VertexOperator jAlphaAlpha;
  jAlphaAlpha.addTerm(new InterProdPartAtVertices(&alphaP));
  decSphere.addMatrixOperator(jAlphaAlpha, 2, 0);

  VertexOperator jHAlphaHAlpha;
  jHAlphaHAlpha.addTerm(new InterProdPartAtVertices(&alphaD));
  decSphere.addMatrixOperator(jHAlphaHAlpha, 2, 1);

  VertexOperator fId;
  fId.addTerm(new IdentityAtVertices(-1.0));
  decSphere.addMatrixOperator(fId, 2, 2);


  decSphere.assembleSystem();
  //cout << decSphere.getSysMat() << endl;
  
  decSphere.solve();

  decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  double errL2Rel_P = decSphere.getEdgeSolution(0).errorL2Rel(alphaP);
  double errMaxRel_P = decSphere.getEdgeSolution(0).errorMaxRel(alphaP);

  double errL2Rel_D = decSphere.getEdgeSolution(1).errorL2Rel(alphaD);
  double errMaxRel_D = decSphere.getEdgeSolution(1).errorMaxRel(alphaD);

  cout << endl;
  cout << "P: RelError L2:  " << errL2Rel_P << endl;
  cout << "P: RelError Max: " << errMaxRel_P << endl;
  cout << endl;
  cout << "D: RelError L2:  " << errL2Rel_D << endl;
  cout << "D: RelError Max: " << errMaxRel_D << endl;


  DofVertexVector magAlpha(edgeMesh, "norm2Alpha");
  magAlpha.interpol(new Norm2DZ());

  double errL2Rel_f = decSphere.getVertexSolution(2).errorL2Rel(magAlpha);
  double errMaxRel_f = decSphere.getVertexSolution(2).errorMaxRel(magAlpha);

  cout << endl;
  cout << "f: RelError L2:  " << errL2Rel_f << endl;
  cout << "f: RelError Max: " << errMaxRel_f << endl;

  AMDiS::finalize();
}
