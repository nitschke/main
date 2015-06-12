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
  F() : AbstractFunction<double, WorldVector<double> >(1) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    //cout << (*timePtr) << endl;
    double pit = M_PI * (*timePtr);
    return (M_PI * cos(pit) + 2.0 * sin(pit)) * x[0];
  }
};

class InitFun : public AbstractFunction<double, WorldVector<double> >,
	  public TimedObject
{
public:
  InitFun() : AbstractFunction<double, WorldVector<double> >(1) {
    srand(12);
  }

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& coords) const 
  {
    //return myrand(0.1,1.5);
    //
    double l = 0.01;
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    double cp4 = cos(M_PI / 4.0);
    double eps = 1.e-16;

    WorldVector<double> n;
    n[0] = x;
    n[1] = y;
    n[2] = z;
    
    WorldVector<double> e;
    e *= 0.0;
    if (abs(y) >= cp4) {
      e[0] = -x;
      e[1] = 0.0;
      e[2] = -z;
    } 
    else if (x >= cp4) {
      e[0] = 0.0;
      e[1] = y;
      e[2] = z;
    } else if (x <= -cp4) {
      e[0] = 0.0;
      e[1] = sin(M_PI * (y - l));
      e[2] = -sin(M_PI *z);
    } else {
      double c = y / cp4;
      e[0] = abs(c) - 1.0;
      e[1] = c;
      e[2] = 0.0; 
    }
    e += -(e*n)*n;
    double norme = sqrt(e*e);
    e *= (norme < eps) ? 0.0 : (1./norme);
    
    WorldVector<double> t;
    t[0] = -y;
    t[1] = x;
    t[2] = 0.0;
    double normt = sqrt(t*t);
    t *= (normt < eps) ? 0.0 : (1./normt);

    double arg = e*t;

    double rval = 0.0;
    if (abs(y) >= cp4) {
      rval = (z < 0) ? (-M_PI + acos(-arg)) : (acos(arg));
    } 
    else if (x >= cp4) {
      rval = (z > 0) ? (-M_PI + acos(-arg)) : (acos(arg));
    } else if (x <= -cp4) {
      rval = (z < 0) ? (-M_PI + acos(-arg)) : (acos(arg));
    } else {
      rval = (false) ? (-M_PI + acos(-arg)) : (acos(arg));
    }

    
    return rval;
  }

private:
  inline double myrand(double a, double b) const{
    return (b-a) * ((double)rand()) / RAND_MAX + a;
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


  InitFun *iFun = new InitFun();
  iFun->setTimePtr(&zero);
  heat.setExactSolution(iFun);

  // ===== start adaption loop =====
  int errorCode = adaptInstat.adapt();

  AMDiS::finalize();



  return errorCode;
}
