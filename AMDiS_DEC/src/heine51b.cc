#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================
class GaussCurv : public AbstractFunction<double, WorldVector<double> >
{
public:
  GaussCurv() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& xx) const 
  {
    double x = xx[0];
    double y = xx[1];
    double z = xx[2];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double tmp = 1.0 - 4.0*(-2.0 + x + x2 + y - 2.0*x*y + y2) * z2;
    return -(1.0 + 2.0*x - 2.0*x2 + 2.0*y - 2.0*y2 + 2.0*(-3.0 + x + y) * z2 ) / (tmp*tmp);
  }
};




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
    rval[2] += 2.0 * x[2]* (2.0 * z2 - x[0] - x[1]);
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
  new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-6, 100);

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

  DOFVector<double> *gaussCurvDV = new DOFVector<double>(sphere.getFeSpace(),"gaussCurv");
  GaussCurv *gc = new GaussCurv();
  WorldVector<double> xtest;
  xtest = 0.0;
  xtest[1] = -1.0;
  cout << (*gc)(xtest) << endl;
  gaussCurvDV->interpol(gc);
  VtkVectorWriter::writeFile(gaussCurvDV, "output/gaussCurv.vtu");

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


