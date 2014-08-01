#include "AMDiS.h"
#include "decOperator.h"
#include "DOFVHelper.h"
#include "torusProjection.h"


using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

class F : public AbstractFunction<double, WorldVector<double> >
{
public:
  F() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& coords) const 
  {
    double r = 0.5;
    double R = 2;
    double x = coords[0];
    double y = coords[2]; //coords change
    double z = coords[1];
    double sx2py2 =  sqrt(x*x + y*y);
    return y / sx2py2 + sx2py2 - R;
  }
};

class GradF : public AbstractFunction< WorldVector<double>, WorldVector<double> >
{
public:
  GradF() : AbstractFunction< WorldVector<double>, WorldVector<double> >(1) {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const
  {
    double r = 0.5;
    double r2 = r*r;
    double R = 2;
    double x = coords[0];
    double y = coords[2]; //coords change y <-> z bzw 2 <-> 1
    double z = coords[1];
    double z2 = z*z;
    double x2py2 =  x*x + y*y;
    double sx2py2 =  sqrt(x2py2);
    WorldVector<double> rvalt;
    rvalt[0] = -x*y;
    rvalt[2] =  x*x;
    rvalt[1] =  0.0;
    rvalt *= 1.0/(sx2py2*x2py2);
    WorldVector<double> rvalp;
    rvalp[0] =  x*z2/sx2py2;
    rvalp[2] =  y*z2/sx2py2;
    rvalp[1] =  - z*(sx2py2 - R);
    rvalp *= 1.0/r2;
    return rvalt + rvalp;
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("torus main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  new TorusProject(1, VOLUME_PROJECTION, 2.0, 0.5);

  // ===== create and init the scalar problem ===== 
  ProblemStat torus("torus");
  torus.initialize(INIT_ALL);



  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("torus->adapt", torus.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("torus->adapt",
					       &torus,
					       adaptInfo);
  
  // Grad(F)
  for (int i = 0; i < 3; i++) {
    // ===== create matrix operator =====
    SimpleDEC *decOperator = new SimpleDEC(torus.getFeSpace(i),torus.getFeSpace(i) );
    torus.addMatrixOperator(decOperator, i, i);

    // ===== create rhs operator =====
    PrimalPrimalGradFunctionDEC *rhsOperator = new PrimalPrimalGradFunctionDEC(i, new F(), torus.getFeSpace(i));
    torus.addVectorOperator(rhsOperator, i);
  }
  // ===== start adaption loop =====
  adapt->adapt();

  DOFVector<WorldVector<double> > solDOFV(torus.getFeSpace(),"solDOFV");
  solDOFV.interpol(new GradF());
  VtkVectorWriter::writeFile(solDOFV, string("output/sol.vtu"));  
  printError(*(torus.getSolution()),0,1,2, solDOFV, "Error");

  DOFVector<double> rhsDofVector(torus.getFeSpace(),"RHS");
  rhsDofVector.interpol(new F());
  VtkVectorWriter::writeFile(rhsDofVector, string("output/torusRHS.vtu"));  


  torus.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


