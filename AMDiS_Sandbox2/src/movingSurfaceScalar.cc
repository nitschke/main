#include "Dec.h"
#include "io/VtkVectorWriter.h"
#include "MeshMover.h"

using namespace std;
using namespace AMDiS;
using namespace dec;


class InitFun : public AbstractFunction<double, WorldVector<double> >
{
public:
  InitFun(double *timePtr) : AbstractFunction<double, WorldVector<double> >(), tptr(timePtr) {
  }

  double operator()(const WorldVector<double>& coords) const 
  {
    double t = *tptr;
    return t*t*coords[0]*coords[0]*coords[2];
  }

private:
  double * tptr;
};

class DtFun : public AbstractFunction<double, WorldVector<double> >
{
public:
  DtFun(double *timePtr) : AbstractFunction<double, WorldVector<double> >(), tptr(timePtr) {
  }

  double operator()(const WorldVector<double>& coords) const 
  {
    double t = *tptr;
    return 2.0 * t * (1.0 + 6.0*t - 8.0*t*t) * coords[0]*coords[0]*coords[2] / (4.0 * t * (1.0 - t) + 1.0);
  }

private:
  double * tptr;
};

class StrechingSphere : public BinaryAbstractFunction<WorldVector<double>, WorldVector<double>, double> {
public:
  
  StrechingSphere() : BinaryAbstractFunction<WorldVector<double>, WorldVector<double>, double>() {}

  WorldVector<double> operator()(const WorldVector<double> &coordRef, const double &t) const {
    WorldVector<double> newCoord = coordRef;
    newCoord[2] *= 4.0 * t * (1.0 - t) + 1.0;
    return newCoord;
  }

};


class MyInstat : public DecProblemInstat {
public:
  MyInstat(DecProblemStat *probStat, const DOFVector<double> &initSol, MeshMover *meshmover)
      : DecProblemInstat(probStat), sol(initSol), mmover(meshmover) {
    rhsFun = new DtFun(getTimePtr());
    rhs = DOFVector<double>(probStat->getMesh()->getFeSpace(), "RHS");
    rhs.interpol(rhsFun);
  
    solFun = new InitFun(getTimePtr());
    sol = DOFVector<double>(probStat->getMesh()->getFeSpace(), "exact Solution");
    sol.interpol(solFun);
  }

  void initTimestep() {
    DecProblemInstat::initTimestep();
    mmover->move(t);
    rhs.interpol(rhsFun);
    sol.interpol(solFun);
  }
  

  void closeTimestep() {
    DecProblemInstat::closeTimestep();
  }

  DOFVector<double>* getSolPtr() {
    return &sol;
  }

  DOFVector<double>* getRHSPtr() {
    return &rhs;
  }

private:
  DOFVector<double> sol;
  DOFVector<double> rhs;
  AbstractFunction<double, WorldVector<double> > *rhsFun;
  AbstractFunction<double, WorldVector<double> > *solFun;
  MeshMover *mmover;
};

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  //WorldVector<double> ballCenter;
  //ballCenter.set(0.0);
  //new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  MeshMover mmover(edgeMesh, new StrechingSphere);

  DOFVector<double> init(sphere.getFeSpace(), "init");
  double zero = 0.0;
  init.interpol(new InitFun(&zero));
  io::VtkVectorWriter::writeFile(init, "output/init.vtu");


  DecProblemStat decSphere(&sphere, edgeMesh);

  MyInstat sphereInstat(&decSphere, init, &mmover);

  VertexOperator Dt;
  Dt.addTerm(new IdentityAtVertices());
  Dt.setUhOld(init);
  decSphere.addMatrixOperator(Dt, 0, 0, sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(Dt, 0, sphereInstat.getInvTauPtr());

  VertexOperator RHS;
  RHS.addTerm(new VertexVecAtVertices(sphereInstat.getRHSPtr()));
  decSphere.addVectorOperator(RHS, 0);

  // exact
  VertexOperator SolLeft;
  SolLeft.addTerm(new IdentityAtVertices());
  decSphere.addMatrixOperator(SolLeft, 1, 1);

  VertexOperator SolRight;
  SolRight.addTerm(new VertexVecAtVertices(sphereInstat.getSolPtr()));
  decSphere.addVectorOperator(SolRight, 1);


  sphereInstat.solve();
}
