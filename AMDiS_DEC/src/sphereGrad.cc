#include "AMDiS.h"
#include "decOperator.h"

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
  
  // Grad(X_2)
  for (int i = 0; i < 3; i++) {
    // ===== create matrix operator =====
    SimpleDEC *decOperator = new SimpleDEC(sphere.getFeSpace(i),sphere.getFeSpace(i) );
    sphere.addMatrixOperator(decOperator, i, i);

    // ===== create rhs operator =====
    PrimalPrimalGradFunctionDEC *rhsOperator = new PrimalPrimalGradFunctionDEC(i, new X(2), sphere.getFeSpace(i));
    sphere.addVectorOperator(rhsOperator, i);
  }
  // Div(X)
  SimpleDEC *decOperator = new SimpleDEC(sphere.getFeSpace(3),sphere.getFeSpace(3) );
  sphere.addMatrixOperator(decOperator, 3, 3);
  for (int i = 0; i < 3; i++) {
    PrimalPrimalGradFunctionDEC *rhsOperator = new PrimalPrimalGradFunctionDEC(i, new X(i), sphere.getFeSpace(3));
    sphere.addVectorOperator(rhsOperator, 3);
  }

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


