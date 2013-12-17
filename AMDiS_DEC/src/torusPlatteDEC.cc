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
    return - (R * (r*R + (r*r + R*R ) * zeta/r) * x[0]) / (r*r*r * sgt*sgt*sgt*sgt);
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("torus main");

  AMDiS::init(argc, argv);

  // ===== no projection, use finalized meshes =====

  // ===== create and init the scalar problem ===== 
  ProblemStat torus("torus");
  torus.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("torus->adapt", torus.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("torus->adapt",
					       &torus,
					       adaptInfo);
  
  // ===== create matrix operator =====
  //Operator laplaceOperator(torus.getFeSpace());
  //laplaceOperator.addTerm(new Simple_SOT(-1.0));
  LBeltramiDEC laplaceOperator_h(torus.getFeSpace(0), torus.getFeSpace(0));
  torus.addMatrixOperator(&laplaceOperator_h, 0, 0);

  LBeltramiDEC laplaceOperator_u(torus.getFeSpace(1), torus.getFeSpace(1));
  torus.addMatrixOperator(&laplaceOperator_u, 1, 1);

  SimpleDEC operator_h(torus.getFeSpace(1), torus.getFeSpace(0));
  operator_h.setFactor(-1.0);
  torus.addMatrixOperator(&operator_h, 1, 0);

  int degree = torus.getFeSpace()->getBasisFcts()->getDegree();

  // ===== create rhs operator =====
  //Operator rhsOperator(torus.getFeSpace());
  //rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  FunctionDEC rhsOperator(torus.getFeSpace(), new F(degree));
  torus.addVectorOperator(&rhsOperator, 0);


  

  // ===== start adaption loop =====
  adapt->adapt();

  //cout << torus.getSystemMatrix(0,0)->getBaseMatrix() << endl;

  DOFVector<double> rhsDofVector(torus.getFeSpace(),"RHS");
  rhsDofVector.interpol(new F(degree));
  VtkVectorWriter::writeFile(rhsDofVector, string("output/torusRHS.vtu"));

  torus.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


