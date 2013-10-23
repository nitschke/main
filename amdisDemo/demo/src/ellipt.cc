#include "AMDiS.h"

#include "boost/date_time/posix_time/posix_time.hpp"

using namespace AMDiS;
using namespace boost::posix_time;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// Dirichlet boundary function
class G : public AbstractFunction<double, WorldVector<double> >
{
public:

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return exp(-10.0 * (x * x));
  }
};

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >
{
public:

  F() : AbstractFunction<double, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    int dow = Global::getGeo(WORLD);
    double r2 = (x * x);
    double ux = exp(-10.0 * r2);
    return -(400.0 * r2 - 20.0 * dow) * ux;
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat ellipt("ellipt");
  ellipt.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo adaptInfo("ellipt->adapt", ellipt.getNumComponents());


  // === create adapt ===
  AdaptStationary adapt("ellipt->adapt", ellipt, adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(ellipt.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  ellipt.addMatrixOperator(matrixOperator, 0, 0);


  // ===== create rhs operator =====
  Operator rhsOperator(ellipt.getFeSpace());
  rhsOperator.addTerm(new CoordsAtQP_ZOT(new F));
  ellipt.addVectorOperator(rhsOperator, 0);
  
  // ===== add boundary conditions =====
  ellipt.addDirichletBC(1, 0, 0, new G);


  // ===== start adaption loop =====
  adapt.adapt();
  
  ellipt.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


