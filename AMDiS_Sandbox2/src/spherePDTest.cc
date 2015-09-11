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
  dz.writeSharpFile("output3/dzSharp.vtu",&sphere);
  DofEdgeVector rotz(edgeMesh, "*dz");
  rotz.setDual(&dz_d);
  rotz.writeSharpFile("output3/rotzSharp.vtu",&sphere);

  DecProblemStat decSphere(&sphere, edgeMesh);

// *** Edge Operators *** //
  double MinusOne = -1.0;

  EdgeOperator LB;
  LB.addTerm(new LaplaceBeltramiAtEdges());
  decSphere.addMatrixOperator(LB,0,0);

  EdgeOperator LCB;
  LCB.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LCB,1,1);

  EdgeOperator Id;
  Id.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(Id,0,0);
  decSphere.addMatrixOperator(Id,1,1);

  EdgeOperator DZ;
  DZ.addTerm(new EdgeVecAtEdges(&dz));
  decSphere.addVectorOperator(DZ, 0);
  decSphere.addVectorOperator(DZ, 1);

  EdgeOperator RotZ;
  RotZ.addTerm(new EdgeVecAtEdges(&rotz));
  decSphere.addVectorOperator(RotZ, 0, &MinusOne);
  decSphere.addVectorOperator(RotZ, 1);

  //EdgeOperator RotZP;
  //RotZP.addTerm(new EdgeVecAtEdges(&rotz));
  //decSphere.addVectorOperator(RotZP, 0, &MinusOne);

  //EdgeOperator RotZD;
  //RotZD.addTerm(new EdgeVecAtEdges(&rotz));
  //decSphere.addVectorOperator(RotZD, 1);


  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();

  dz.writeSharpFile("output3/dzSharp.vtu",&sphere);
  rotz.writeSharpFile("output3/rotzSharp.vtu",&sphere);
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


