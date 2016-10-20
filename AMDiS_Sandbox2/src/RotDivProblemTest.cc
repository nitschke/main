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

class f : public AbstractFunction<double, WorldVector<double> > {
  public:
  f() : AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& v) const {
    return v[2];
  }
};

class LBf_Ellipsoid : public AbstractFunction<double, WorldVector<double> > {
  public:
  LBf_Ellipsoid() : AbstractFunction<double, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& v) const {
    double y2 = v[1]*v[1];
    double z = v[2];
    double z2 = z*z;
    double nenSqr = 81. + 972.*y2 - 20.*z2;
    return 288.*z*(-45. - 54.*y2 +10.*z2)/(nenSqr*nenSqr);
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

// Definition of p //
  
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dff = new Df();
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *rotff = new RotZ_Sphere();

  DofEdgeVector dfP(edgeMesh, "df");
  dfP.set(dff);
  DofEdgeVector dfD(edgeMesh, "Rotf");
  dfD.interpolGL4(rotff, sproj.getProjection(), sproj.getJProjection());

  DofEdgeVector p(edgeMesh, "p");
  p = dfP + dfD;

// Definition of f // 
  

  AbstractFunction<double, WorldVector<double> > *ff = new f();
  DofVertexVector fvec(edgeMesh, "f");
  fvec.interpol(ff);


  // +p (= *dz + dz)
  EdgeOperator pOp;
  pOp.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(pOp, 0, 0);

  // + *df (= *dz)
  EdgeOperator rotfOp; 
  rotfOp.addTerm(new RotAtEdges(1.0));
  decSphere.addMatrixOperator(rotfOp, 0, 1);

  // = q (= 2*dz  = 2p)
  EdgeOperator qOp;
  qOp.addTerm(new EdgeVecAtEdges(&dfD, 2.0));
  qOp.addTerm(new EdgeVecAtEdges(&dfP, 1.0));
  decSphere.addVectorOperator(qOp, 0);

  // div(p) (= *d*(*dp + dp)= -*ddp + LB(z) = LB(z)) 
  //VertexOperator divpOp;
  //divpOp.addTerm(new DivAtVertices());
  //decSphere.addMatrixOperator(divpOp, 1, 0);
  VertexOperator divpOp;
  divpOp.addTerm(new VertexVecAtVertices(&fvec, 2.0));
  decSphere.addVectorOperator(divpOp, 1);

  // f ( = z)
  VertexOperator fOp;
  fOp.addTerm(new IdentityAtVertices());
  decSphere.addMatrixOperator(fOp, 1, 1);

  // = h ( = LB(z) + z= z = f)
  DofVertexVector lbfvec(edgeMesh, "lbf");
  //lbfvec.interpol(new LBf_Ellipsoid());
  lbfvec = - 2.0 * fvec; //on Sphere
  
  VertexOperator hOp;
  hOp.addTerm(new VertexVecAtVertices(&fvec, 1.0));
  hOp.addTerm(new VertexVecAtVertices(&lbfvec, 1.0));
  decSphere.addVectorOperator(hOp, 1);


  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  double errL2Rel_p = decSphere.getEdgeSolution(0).errorL2Rel(p);
  double errMaxRel_p = decSphere.getEdgeSolution(0).errorMaxRel(p);

  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  cout << endl;
  cout << "p: RelError L2:  " << errL2Rel_p << endl;
  cout << "p: RelError Max: " << errMaxRel_p << endl;

  DofVertexVector fhvec = decSphere.getVertexSolution(1);
  double errMaxRel_f = fhvec.errorMaxRel(fvec);
  double errL2Rel_f = fhvec.errorL2Rel(fvec);
  cout << endl;
  cout << "f: RelError L2:  " << errL2Rel_f << endl;
  cout << "f: RelError Max: " << errMaxRel_f << endl;

  DofVertexVector *lbhfvec = fvec.laplace();
  double errMaxRel_lbf = lbhfvec->errorMaxRel(lbfvec);
  double errL2Rel_lbf = lbhfvec->errorL2Rel(lbfvec);
  cout << endl;
  cout << "lbf: RelError L2:  " << errL2Rel_lbf << endl;
  cout << "lbf: RelError Max: " << errMaxRel_lbf << endl;

  //cout << endl;
  //cout << "maxDiameter,maxLength,errorL2Rel,errorMaxRel" << endl;
  //cout << hdia << "," << hlen << "," << errL2Rel << ", " << errMaxRel << endl;

  AMDiS::finalize();
}
