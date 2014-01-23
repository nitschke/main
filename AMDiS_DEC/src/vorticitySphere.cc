#include "AMDiS.h"
#include "decOperator.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >,
	  public TimedObject
{
public:
  F(double viscosity) : AbstractFunction<double, WorldVector<double> >(1), nu(viscosity) {}

  double operator()(const WorldVector<double>& x) const 
  {
    double g = 1.0 - 2.0 * (*timePtr);
    double dg = -2.0;
    double x2 = x[0] * x[0];
    double y2 = x[1] * x[1];
    return -112.0 * x[0] * x[1] * x[2] * (x2 - 3.0 * y2) * (3.0 * x2 - y2) * (dg + 54.0 * nu * g);
  }

protected:
  double nu

};

class Psi : public AbstractFunction<double, WorldVector<double> >,
	  public TimedObject
{
public:
  Psi() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& x) const 
  {
    double g = 1.0 - 2.0 * (*timePtr);
    double x2 = x[0] * x[0];
    double y2 = x[1] * x[1];
    return -2.0 * x[0] * x[1] * x[2] * (x2 - 3.0 * y2) * (3.0 * x2 - y2) * dg;
  }
};


class Phi : public AbstractFunction<double, WorldVector<double> >,
	  public TimedObject
{
public:
  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& x) const 
  {
    double g = 1.0 - 2.0 * (*timePtr);
    double x2 = x[0] * x[0];
    double y2 = x[1] * x[1];
    return -112.0 * x[0] * x[1] * x[2] * (x2 - 3.0 * y2) * (3.0 * x2 - y2) * dg;
  }
};


// ===========================================================================
// ===== instationary problem ================================================
// ===========================================================================

/// Instationary problem
class Heat : public ProblemInstat
{
public:
  Heat(ProblemStat &heatSpace)
    : ProblemInstat("heat", heatSpace)
  {
    // init theta scheme
    theta = -1.0;
    Parameters::get(name + "->theta", theta);
    TEST_EXIT(theta >= 0)("theta not set!\n");
    if (theta == 0.0) {
      WARNING("You are using the explicit Euler scheme\n");
      WARNING("Use a sufficiently small time step size!!!\n");
    }
    MSG("theta = %f\n", theta);
    theta1 = theta - 1.0;
  }

  // ===== ProblemInstatBase methods ===================================

  /// set the time in all needed functions!
  void setTime(AdaptInfo *adaptInfo) 
  {
    ProblemInstat::setTime(adaptInfo);
    rhsTime = adaptInfo->getTime() - (1 - theta) * adaptInfo->getTimestep();
  }

  void closeTimestep(AdaptInfo *adaptInfo) 
  {
    ProblemInstat::closeTimestep(adaptInfo);
    WAIT;
    //cout << problemStat->getSystemMatrix(0,0)->getBaseMatrix() << endl;
  }

  // ===== initial problem methods =====================================

  /// Used by \ref problemInitial to solve the system of the initial problem
  void solve(AdaptInfo *adaptInfo) 
  {
    problemStat->getSolution(0)->interpol(exactSolution);
  }

  /// Used by \ref problemInitial to do error estimation for the initial problem.
  void estimate(AdaptInfo *adaptInfo) 
  {
    double errMax, errSum;

    errSum = Error<double>::L2Err(*exactSolution,
				  *(problemStat->getSolution(0)), 0, &errMax, false);
    adaptInfo->setEstSum(errSum, 0);
    adaptInfo->setEstMax(errMax, 0);
  }

  // ===== setting methods ===============================================

  /// Sets \ref exactSolution;
  void setExactSolution(AbstractFunction<double, WorldVector<double> > *fct) 
  {
    exactSolution = fct;
  } 

  // ===== getting methods ===============================================

  /// Returns pointer to \ref theta.
  double *getThetaPtr() 
  { 
    return &theta; 
  }

  /// Returns pointer to \ref theta1.
  double *getTheta1Ptr() 
  { 
    return &theta1; 
  }

  /// Returns pointer to \ref rhsTime.
  double *getRhsTimePtr() 
  { 
    return &rhsTime; 
  }

private:
  /// Used for theta scheme.
  double theta;

  /// theta - 1
  double theta1;

  /// time for right hand side functions.
  double rhsTime;

  /// Pointer to boundary function. Needed for initial problem.
  AbstractFunction<double, WorldVector<double> > *exactSolution;
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char** argv)
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init stationary problem =====
  ProblemStat heatSpace("heat->space");
  heatSpace.initialize(INIT_ALL);


  // ===== create instationary problem =====
  Heat heat(heatSpace);
  heat.initialize(INIT_ALL);

  // create adapt info
  AdaptInfo adaptInfo("heat->adapt", heatSpace.getNumComponents());

  // create initial adapt info
  AdaptInfo adaptInfoInitial("heat->initial->adapt");

  // create instationary adapt
  AdaptInstationary adaptInstat("heat->adapt",
				heatSpace,
				adaptInfo,
				heat,
				adaptInfoInitial);


  // ===== create rhs functions =====
  F *rhsFct = new F();
  rhsFct->setTimePtr(heat.getRhsTimePtr());


  // ===== create operators =====
  double one = 1.0;
  double zero = 0.0;

  // create laplace
  //Operator A(heatSpace.getFeSpace());
  //A.addTerm(new Simple_SOT);
  LBeltramiDEC A(heatSpace.getFeSpace());
  A.setFactor(-1.0);
  //A.setUhOld(heat.getOldSolution(0));
  //if (*(heat.getThetaPtr()) != 0.0)
    heatSpace.addMatrixOperator(A, 0, 0, heat.getThetaPtr());
  //if (*(heat.getTheta1Ptr()) != 0.0)
  //  heatSpace.addVectorOperator(A, 0, heat.getTheta1Ptr(), &zero);

  // create zero order operator
  //Operator C(heatSpace.getFeSpace());
  //C.addTerm(new Simple_ZOT);
  double *invTau = heat.getInvTau();
  //invTau = &zero;
  SimpleDEC C(heatSpace.getFeSpace());
  C.setUhOld(heat.getOldSolution(0));
  heatSpace.addMatrixOperator(C, 0, 0, invTau);
  heatSpace.addVectorOperator(C, 0, invTau);

  // create RHS operator
  //Operator F(heatSpace.getFeSpace());
  //F.addTerm(new CoordsAtQP_ZOT(rhsFct));
  FunctionDEC F(rhsFct, heatSpace.getFeSpace());
  heatSpace.addVectorOperator(F, 0);

  InitFun *iFun = new InitFun();
  iFun->setTimePtr(&zero);
  heat.setExactSolution(iFun);

  // ===== start adaption loop =====
  int errorCode = adaptInstat.adapt();

  AMDiS::finalize();



  return errorCode;
}
