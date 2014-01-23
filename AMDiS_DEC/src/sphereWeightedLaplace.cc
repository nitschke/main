#include "AMDiS.h"
#include "decOperator.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >
{
public:
  F(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return -2.0 * x[0];
    //return -6.0 * x[0] * x[2];
  }
};

class W : public AbstractFunction<double, WorldVector<double> >
{
public:
  W(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return 1.0 + 10.0 * x[1]*x[1];
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
  sphere.setWriteAsmInfo(true);
  sphere.setAssembleMatrixOnlyOnce(0, 0, false);



  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  
  int degree = 1;

  LBeltramiWeightedDEC decOperator(new W(degree), sphere.getFeSpace());
  sphere.addMatrixOperator(&decOperator, 0, 0);

  FunctionDEC rhsOperator(new F(degree), sphere.getFeSpace());
  sphere.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  //cout << sphere.getSystemMatrix(0,0)->getBaseMatrix() << endl;
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;
  //sphere.getSystemMatrix(0,0)->calculateNnz();
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


