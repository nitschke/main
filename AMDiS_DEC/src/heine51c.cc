#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"
#include "meshCorrector.h"
#include "MeshHelper.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

class X : public AbstractFunction<double, WorldVector<double> >
{
public:
  X(int i_) : AbstractFunction<double, WorldVector<double> >(1), i(i_) {}

  double operator()(const WorldVector<double>& x) const 
  {
    return x[i];
  }

protected:
  int i;
};

class Phi : public AbstractFunction<double, WorldVector<double> >
{
public:
  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& x) const 
  {
    return 0.5 * (x[0]*x[0] + 4.0*x[1]*x[1] + (4.0/9.0)*x[2]*x[2] - 1.0);
  }
};

class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}

  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> rval(x);
    rval[1] *= 4.0;
    rval[2] *= 4.0/9.0;
    return rval;
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
  new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-6);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);



  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  SimpleDEC *gaussCurv = new SimpleDEC(sphere.getFeSpace(0),sphere.getFeSpace(0) );
  sphere.addMatrixOperator(gaussCurv, 0, 0);

  GaussCurvatureDEC *rhs = new GaussCurvatureDEC(sphere.getFeSpace(0));
  sphere.addVectorOperator(rhs,0);

  // ===== start adaption loop =====
  adapt->adapt();

  DOFVector<WorldVector<double> > forces = getConnectionForces(sphere.getFeSpace());
  VtkVectorWriter::writeFile(forces, string("output/ConForces.vtu"));

  MeshCorrector mc(sphere.getFeSpace());
  int n = 10000;
  for (int i = 0; i < n; i++) {
    mc.oneIteration(1.0);
  }
  sphere.setFeSpace(mc.getFeSpace());
  DOFVector<double> phiNew(mc.getFeSpace(), "phiNew");
  phiNew.interpol(new Phi());
  VtkVectorWriter::writeFile(phiNew, string("output/phiNew.vtu"));

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


