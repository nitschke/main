#include "Dec.h"
#include "io/VtkVectorWriter.h"
#include "torusProjection.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

// <dxyz,[p,q]>
class DXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
  }
};

class RotXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    WorldVector<double> rot;
    double x = coords[0];
    double y = coords[2]; //coords change 2 <-> 1
    double z = coords[1];
    rot[0] = 0.015625*x*pow(pow(x,2.) + pow(y,2.),-1.5)*(8.*pow(x,4.)*
       (-2. + pow(pow(x,2.) + pow(y,2.),0.5) + pow(y,2.)*(10. + 3.*pow(pow(x,2.) + pow(y,2.),0.5))) + 
      2.*pow(x,2.)*(30. + pow(z,2.)*(8. - 12.*pow(pow(x,2.) + pow(y,2.),0.5)) - 17.*pow(pow(x,2.) + pow(y,2.),0.5) + 
         8.*pow(y,4.)*(10. + 3.*pow(pow(x,2.) + pow(y,2.),0.5)) - 
         5.*pow(y,2.)*(34. + 5.*pow(pow(x,2.) + pow(y,2.),0.5) + 24.*pow(z,2.)*(2. + pow(pow(x,2.) + pow(y,2.),0.5)))) + 
      pow(y,2.)*(90. - 51.*pow(pow(x,2.) + pow(y,2.),0.5) + 40.*pow(z,4.)*(2. + 3.*pow(pow(x,2.) + pow(y,2.),0.5)) + 
         8.*pow(y,4.)*(10. + 3.*pow(pow(x,2.) + pow(y,2.),0.5)) + 6.*pow(z,2.)*(54. + 29.*pow(pow(x,2.) + pow(y,2.),0.5)) - 
         2.*pow(y,2.)*(162. + 29.*pow(pow(x,2.) + pow(y,2.),0.5) + 120.*pow(z,2.)*(2. + pow(pow(x,2.) + pow(y,2.),0.5)))));
    rot[1] = 0.03125*(x - 1.*y)*(x + y)*z*pow(pow(x,2.) + pow(y,2.),-1.)*
    (-15. + 12.*pow(x,2.) + 12.*pow(y,2.) - 4.*pow(z,2.) - 16.*pow(pow(x,2.) + pow(y,2.),0.5));
    rot[2] = -0.015625*y*pow(pow(x,2.) + pow(y,2.),-1.5)*(2.*pow(y,2.)*
       (30. + pow(z,2.)*(8. - 12.*pow(pow(x,2.) + pow(y,2.),0.5)) + 4.*pow(y,2.)*(-2. + pow(pow(x,2.) + pow(y,2.),0.5)) - 
         17.*pow(pow(x,2.) + pow(y,2.),0.5)) + 8.*pow(x,6.)*(10. + 3.*pow(pow(x,2.) + pow(y,2.),0.5)) + 
      2.*pow(x,4.)*(-162. - 29.*pow(pow(x,2.) + pow(y,2.),0.5) - 120.*pow(z,2.)*(2. + pow(pow(x,2.) + pow(y,2.),0.5)) + 
         8.*pow(y,2.)*(10. + 3.*pow(pow(x,2.) + pow(y,2.),0.5))) + 
      pow(x,2.)*(90. - 51.*pow(pow(x,2.) + pow(y,2.),0.5) + 40.*pow(z,4.)*(2. + 3.*pow(pow(x,2.) + pow(y,2.),0.5)) + 
         8.*pow(y,4.)*(10. + 3.*pow(pow(x,2.) + pow(y,2.),0.5)) + 6.*pow(z,2.)*(54. + 29.*pow(pow(x,2.) + pow(y,2.),0.5)) - 
         10.*pow(y,2.)*(34. + 5.*pow(pow(x,2.) + pow(y,2.),0.5) + 24.*pow(z,2.)*(2. + pow(pow(x,2.) + pow(y,2.),0.5)))));
    return  rot * vec;
  }
};

// <dt,[p,q]>
class DT_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DT_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double tp = (p[2] < 0) ? (M_PI + acos(-p[0]/sqrt(p[0]*p[0]+p[2]*p[2]))) : (acos(p[0]/sqrt(p[0]*p[0]+p[2]*p[2])));
    double tq = (q[2] < 0) ? (M_PI + acos(-q[0]/sqrt(q[0]*q[0]+q[2]*q[2]))) : (acos(q[0]/sqrt(q[0]*q[0]+q[2]*q[2])));
    double rval = tq - tp;
    if (abs(rval) > M_PI) rval += (tp > M_PI) ? (2.0*M_PI) : (-2.0*M_PI);
    return  rval;
  }
};


// e.g. Laplace-Beltrami of exact forms or Laplace-CoBeltrami of co-exact forms
class Zero : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Zero() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    return 0.0;
  }
};


////projection
//class Proj : public AbstractFunction<WorldVector<double>, WorldVector<double> > {
//
//  public:
//  Proj() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}
//
//  /// Implementation of AbstractFunction::operator().
//  WorldVector<double> operator()(const WorldVector<double>& x) const 
//  {
//  }
//};
//
////jacobi matrix of projection
//class JProj : public AbstractFunction<WorldMatrix<double>, WorldVector<double> > {
//
//  public:
//  JProj() : AbstractFunction<WorldMatrix<double>, WorldVector<double> >() {}
//
//  /// Implementation of AbstractFunction::operator().
//  WorldMatrix<double> operator()(const WorldVector<double>& x) const 
//  {
//  }
//
//
//};




// ===========================================================================
// ===== main program ========================================================
// ===========================================================================


// harmonischer anteil bleibt in allen loesungen erhalten
int main(int argc, char* argv[])
{
  FUNCNAME("torus main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  new TorusProject(1, VOLUME_PROJECTION, 2.0, 0.5);

  // ===== create and init the scalar problem ===== 
  ProblemStat torus("torus");
  torus.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(torus.getFeSpace());
  
  // div-free
  DofEdgeVector rotxyz(edgeMesh, "initDivFree");
  rotxyz.interpolLinTrapz(new RotXYZ());
  DOFVector< WorldVector<double> > rotxyzSharp = rotxyz.getSharpFaceAverage();
  io::VtkVectorWriter::writeFile(rotxyzSharp, "output/initDivFree.vtu");


  // rot-free
  DofEdgeVector dxyz(edgeMesh, "initRotFree");
  dxyz.set(new DXYZ_d());
  DOFVector< WorldVector<double> > dxyzSharp = dxyz.getSharpFaceAverage();
  io::VtkVectorWriter::writeFile(dxyzSharp, "output/initRotFree.vtu");

  // harmonic
  DofEdgeVector dt(edgeMesh, "initHarm");
  dt.set(new DT_d());
  //dt.writeFile("output/initHarm.vtu");
  DOFVector< WorldVector<double> > dtSharp = dt.getSharpFaceAverage();
  io::VtkVectorWriter::writeFile(dtSharp, "output/initHarmSharp.vtu");




  DofEdgeVector initSol = dxyz + rotxyz + 1.0*dt;
  initSol.setName("initSol");
  initSol.writeFile("output/initSol.vtu");
  DOFVector< WorldVector<double> > initSolSharp = initSol.getSharpFaceAverage();
  io::VtkVectorWriter::writeFile(initSolSharp, "output/initSolSharp.vtu");


  DecProblemStat dectorus(&torus, edgeMesh);

  DecProblemInstat torusInstat(&dectorus);

  double minusOne = -1.0;
  // -Beltrami
  EdgeOperator Beltrami;
  Beltrami.addTerm(new LaplaceBeltramiAtEdges());
  dectorus.addMatrixOperator(Beltrami, 0, 0, &minusOne);
  // -CoBeltrami
  EdgeOperator CoBeltrami;
  CoBeltrami.addTerm(new LaplaceCoBeltramiAtEdges());
  dectorus.addMatrixOperator(CoBeltrami, 1, 1, &minusOne);
  //// DeRham
  //dectorus.addMatrixOperator(Beltrami, 2, 2, &minusOne);
  //dectorus.addMatrixOperator(CoBeltrami, 2, 2, &minusOne);

  // time derivatives approx
  EdgeOperator I0Operator;
  I0Operator.addTerm(new IdentityAtEdges());
  I0Operator.setUhOld(initSol);
  dectorus.addMatrixOperator(I0Operator, 0, 0,torusInstat.getInvTauPtr());
  dectorus.addVectorOperator(I0Operator, 0,torusInstat.getInvTauPtr());

  EdgeOperator I1Operator;
  I1Operator.addTerm(new IdentityAtEdges());
  I1Operator.setUhOld(initSol);
  dectorus.addMatrixOperator(I1Operator, 1, 1,torusInstat.getInvTauPtr());
  dectorus.addVectorOperator(I1Operator, 1,torusInstat.getInvTauPtr());

  //EdgeOperator I2Operator;
  //I2Operator.addTerm(new IdentityAtEdges());
  //I2Operator.setUhOld(initSol);
  //dectorus.addMatrixOperator(I2Operator, 2, 2,torusInstat.getInvTauPtr());
  //dectorus.addVectorOperator(I2Operator, 2,torusInstat.getInvTauPtr());


  EdgeOperator Id;
  Id.addTerm(new IdentityAtEdges());
  dectorus.addMatrixOperator(Id, 2, 0);
  dectorus.addMatrixOperator(Id, 2, 1);
  dectorus.addMatrixOperator(Id, 2, 2);

  EdgeOperator Alpha0;
  Alpha0.addTerm(new EdgeVecAtEdges(&initSol));
  dectorus.addVectorOperator(Alpha0,2);


  torusInstat.solve();
  



  //// === create adapt info ===
  //AdaptInfo *adaptInfo = new AdaptInfo("torus->adapt", torus.getNumComponents());

  //// === create adapt ===
  //AdaptStationary *adapt = new AdaptStationary("torus->adapt",
	//				       &torus,
	//				       adaptInfo);
  

  // ===== start adaption loop =====
  //adapt->adapt();

  //torus.writeFiles(adaptInfo, true);
  
  AMDiS::finalize();
}


