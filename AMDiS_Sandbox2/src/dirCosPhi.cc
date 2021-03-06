#include "Dec.h"
#include "io/VtkVectorWriter.h"

using namespace std;
using namespace AMDiS;
using namespace dec;


class InitFun : public AbstractFunction<double, WorldVector<double> >
{
public:
  InitFun() : AbstractFunction<double, WorldVector<double> >() {
    srand(12);
  }

  double operator()(const WorldVector<double>& coords) const 
  {
    return myrand(-0.9,0.9);
  }

private:
  inline double myrand(double a, double b) const{
    return (b-a) * ((double)rand()) / RAND_MAX + a;
  }
};

class LBFun : public AbstractFunction<double, double> {
public:
  LBFun() :  AbstractFunction<double, double >() {}

  double operator()(const double &c) const {
    return 1.0 + c / (1.0 - c*c);
  }
};


class MyInstat : public DecProblemInstat {
public:
  MyInstat(DecProblemStat *probStat, const DOFVector<double> &initSol)
      : DecProblemInstat(probStat), sol(initSol) {};

  void closeTimestep() {
    DecProblemInstat::closeTimestep();
    sol = statProb->getVertexSolution(0);
  }

  DOFVector<double>* getSolPtr() {
    return &sol;
  }

private:
  DOFVector<double> sol;
};

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

  DOFVector<double> init(sphere.getFeSpace(), "init");
  init.interpol(new InitFun());
  io::VtkVectorWriter::writeFile(init, "output/init.vtu");


  DecProblemStat decSphere(&sphere, edgeMesh);

  MyInstat sphereInstat(&decSphere, init);

  VertexOperator Laplace;
  Laplace.addTerm(new VertexVecLaplaceBeltramiAtVertices(sphereInstat.getSolPtr(), new LBFun(), -1.0));
  decSphere.addMatrixOperator(Laplace, 0, 0);

  VertexOperator Dt;
  Dt.addTerm(new IdentityAtVertices());
  Dt.setUhOld(init);
  decSphere.addMatrixOperator(Dt, 0, 0, sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(Dt, 0, sphereInstat.getInvTauPtr());

  sphereInstat.solve();
}
