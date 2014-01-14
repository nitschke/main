#include "AMDiS.h"
#include "decOperator.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

class Phi : public AbstractFunction<double, WorldVector<double> >
{
public:
  Phi(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return x[0] * x[2];
    //return 1.0;
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



  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  // ===== create matrix operator =====
  JacobianDEC decOperator(new Phi(1), sphere.getFeSpace());
  sphere.addMatrixOperator(&decOperator, 0, 0);

  // ===== create rhs operator =====
  //FunctionDEC rhsOperator(new F(1), sphere.getFeSpace());
  //sphere.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  cout << "SysMat:" << endl;
  cout << sphere.getSystemMatrix(0,0)->getBaseMatrix() << endl;
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;
  //sphere.getSystemMatrix(0,0)->calculateNnz();
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


