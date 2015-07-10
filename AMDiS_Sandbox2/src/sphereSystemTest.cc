#include "Dec.h"
#include "io/VtkVectorWriter.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

// <*dz,[p,q]> correlate to vec [y, -x, 0]
class Rotz_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Rotz_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return -2.*(p2*q1 - 1.*p1*q2)*atan((-1. + p1*q1 + p2*q2 + p3*q3)*pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + 
       pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5))*
   pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + 
     pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5);
  }
};



// <dz,[p,q]> coorelate to vec [-xz, -yz, 1-z^2]
class DZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[2] - p[2];
  }
};

// <df,[p,q]> with f = z * (3 - z^2/3) -> LdRham(dz) + ||dz||^2 * dz 
class DF_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DF_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  f(q[2]) - f(p[2]);
  }

private:
   inline double f(const double &z) const {
    return z * (3.0 - z * z / 3.0);
   }

};


// <df,[p,q]> with f = z - eps*sqrt(1-z^2) 
class DZDisturbed_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZDisturbed_d(double eps_) : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >(), eps(eps_) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  f(q[2]) - f(p[2]);
  }

private:
   inline double f(const double &z) const{
    return z - eps*sqrt(1.0 - z*z);
   }

   double eps;

};

class Norm2 : public TertiaryAbstractFunction<double,double,double, EdgeElement> {
public:
  Norm2() : TertiaryAbstractFunction<double,double,double, EdgeElement>() {}

  double operator()(const double &primal, const double &dual, const EdgeElement &eel) const {
    double lenP = eel.infoLeft->getEdgeLen(eel.dofEdge);
    return (primal*primal + dual*dual) / (lenP * lenP);
  }
};




// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());
  
  DZ_d dz_d;
  DofEdgeVector dz(edgeMesh, "dz");
  dz.set(&dz_d);
  dz.writeSharpFile("output/dzSharp.vtu",&sphere);
  DofEdgeVector rotz(edgeMesh, "*dz");
  rotz.setDual(&dz_d);
  rotz.writeSharpFile("output/rotzSharp.vtu",&sphere);

  DZDisturbed_d dzd_d(0.1);
  DofEdgeVector dzd(edgeMesh, "dz Disturbed");
  dzd.set(&dzd_d);
  dzd.writeSharpFile("output/dzDisturbedSharp.vtu",&sphere);
  DofEdgeVector rotzd(edgeMesh, "*dz Disturbed");
  rotzd.setDual(&dzd_d);
  rotzd.writeSharpFile("output/rotzDisturbedSharp.vtu",&sphere);

  DecProblemStat decSphere(&sphere, edgeMesh);

// *** Edge Operators *** //
  double MinusOne = -1.0;
 // LaplaceDeRham
  EdgeOperator LaplaceOperator;
  LaplaceOperator.addTerm(new LaplaceBeltramiAtEdges());
  LaplaceOperator.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceOperator, 0, 0, &MinusOne);
  decSphere.addMatrixOperator(LaplaceOperator, 1, 1, &MinusOne);

  EdgeOperator IPOperator;
  IPOperator.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(IPOperator, 0, 1);
  EdgeOperator IDOperator;
  IDOperator.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(IDOperator, 1, 0, &MinusOne);

  EdgeOperator RHSPrimalOperator; // 2dz - *dz
  RHSPrimalOperator.addTerm(new EdgeVecAtEdges(&dz, 2.0));
  RHSPrimalOperator.addTerm(new EdgeVecAtEdges(&rotz, 1.0));
  decSphere.addVectorOperator(RHSPrimalOperator, 0);

  EdgeOperator RHSDualOperator; // -dz + 2(*dz)
  RHSDualOperator.addTerm(new EdgeVecAtEdges(&dz, -1.0));
  RHSDualOperator.addTerm(new EdgeVecAtEdges(&rotz, 2.0));
  decSphere.addVectorOperator(RHSDualOperator, 1);



  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();

  //using namespace mtl;
  //dense2D<double> mat(60,60);
  //mat = decSphere.getSysMat();
  //edgeMesh->printVolInfos(true,true);
  //cout << sub_matrix(mat,0,30,0,30) << endl << endl;
  //cout << sub_matrix(mat,0,30,30,60) << endl << endl;
  //cout << sub_matrix(mat,30,60,0,30) << endl << endl;
  //cout << sub_matrix(mat,30,60,30,60) << endl << endl;
  //cout << decSphere.getRhs() << endl;



  
  AMDiS::finalize();
}


