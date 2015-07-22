#include "Dec.h"
#include "io/VtkVectorWriter.h"
#include "MeshMover.h"

using namespace std;
using namespace AMDiS;
using namespace dec;


class InitFun : public AbstractFunction<double, WorldVector<double> >
{
public:
  InitFun() : AbstractFunction<double, WorldVector<double> >() {
  }

  double operator()(const WorldVector<double>& coords) const 
  {
    return coords[0]*coords[0]*coords[2];
  }
};

class DtFun : public AbstractFunction<double, WorldVector<double> >
{
public:
  DtFun(double *timePtr) : AbstractFunction<double, WorldVector<double> >(), tptr(timePtr) {
  }

  double operator()(const WorldVector<double>& coords) const 
  {
    double t = *tptr;
    return 4.0 * (1.0 - 2.0 * t) * coords[0]*coords[0]*coords[2] / (4.0 * t * (1.0 - t) + 1.0);
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
  }

  void initTimestep() {
    DecProblemInstat::initTimestep();
    mmover->move(t);
    rhs.interpol(rhsFun);
  }
  

  void closeTimestep() {
    DecProblemInstat::closeTimestep();
    sol = statProb->getVertexSolution(0);
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
  init.interpol(new InitFun());
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

  sphereInstat.solve();
}
