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
    double r = 0.5;
    double R = 2.0;
    double zeta =  (abs(x[1]) < r) ? (sqrt(r*r - x[1]*x[1])) : 0.0;
    if (x[0]*x[0] + x[2]*x[2] < R*R) zeta *= -1.0;
    double sgt = R + zeta;
    return 4.0 * x[0];
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== no projection, use finalized meshes =====

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
  //Operator laplaceOperator(sphere.getFeSpace());
  //laplaceOperator.addTerm(new Simple_SOT(-1.0));
  LBeltramiDEC laplaceOperator_h(sphere.getFeSpace(0), sphere.getFeSpace(0));
  sphere.addMatrixOperator(&laplaceOperator_h, 0, 0);

  LBeltramiDEC laplaceOperator_u(sphere.getFeSpace(1), sphere.getFeSpace(1));
  sphere.addMatrixOperator(&laplaceOperator_u, 1, 1);

  SimpleDEC operator_h(sphere.getFeSpace(1), sphere.getFeSpace(0));
  operator_h.setFactor(-1.0);
  sphere.addMatrixOperator(&operator_h, 1, 0);

  int degree = sphere.getFeSpace()->getBasisFcts()->getDegree();

  // ===== create rhs operator =====
  //Operator rhsOperator(sphere.getFeSpace());
  //rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  FunctionDEC rhsOperator(new F(degree), sphere.getFeSpace());
  sphere.addVectorOperator(&rhsOperator, 0);


  

  // ===== start adaption loop =====
  adapt->adapt();

  //cout << sphere.getSystemMatrix(0,0)->getBaseMatrix() << endl;

  DOFVector<double> rhsDofVector(sphere.getFeSpace(),"RHS");
  rhsDofVector.interpol(new F(degree));
  VtkVectorWriter::writeFile(rhsDofVector, string("output/sphereRHS.vtu"));

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


