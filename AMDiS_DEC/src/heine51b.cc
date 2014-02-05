#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"

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
    double z2 = x[2] * x[2];
    double xMz2 = x[0] - z2;
    double yMz2 = x[1] - z2;
    return 0.5 * (xMz2 * xMz2 + yMz2 * yMz2 + z2 - 1.0);
  }
};

class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}

  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> rval(x);
    double z2 = x[2] * x[2];
    rval[0] -= z2;
    rval[1] -= z2;
    rval[2] += 2.0 * x[2]* (z2 - x[0] - x[1]);
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
  new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-2);

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


  //sphere.getRhsVector()->print();

  //cout << sphere.getSystemMatrix(0,0)->getBaseMatrix() << endl;
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;
  //sphere.getSystemMatrix(0,0)->calculateNnz();
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


