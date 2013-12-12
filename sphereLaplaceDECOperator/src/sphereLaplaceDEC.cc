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
    return (x[0] > 0.5) ? 1.0 : 0.0;
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



  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  // ===== create matrix operator =====
  //Operator matrixOperator(sphere.getFeSpace());
  //matrixOperator.addTerm(new Simple_SOT);
  //sphere.addMatrixOperator(&matrixOperator, 0, 0);
  LBeltramiDEC decOperator(sphere.getFeSpace());
  sphere.addMatrixOperator(&decOperator, 0, 0);

  int degree = sphere.getFeSpace()->getBasisFcts()->getDegree();

  // ===== create rhs operator =====
  //Operator rhsOperator(sphere.getFeSpace());
  FunctionDEC rhsOperator(sphere.getFeSpace(), new F(degree));


  //rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  sphere.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  //cout << sphere.getSystemMatrix(0,0)->getBaseMatrix() << endl;

  
  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


