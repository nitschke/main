#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"
#include "meshCorrector.h"
#include "MeshHelper.h"
#include "WorldVectorHelper.h"
#include "MatrixHelper.h"
#include "DOFVHelper.h"

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

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  for (int i = 0; i < 3; i++) {
    // ===== create matrix operator =====
    sphere.addMatrixOperator(new SimpleDEC(sphere.getFeSpace(i), sphere.getFeSpace(i)), i, i);
    // ===== create rhs operator =====
    sphere.addVectorOperator(new LBeltramiInteriorFunctionDEC(new X(i), sphere.getFeSpace(i)), i);
  }
  
  // GaussCurv: geometric operator
  SimpleDEC *gaussCurv = new SimpleDEC(sphere.getFeSpace(3),sphere.getFeSpace(3) );
  sphere.addMatrixOperator(gaussCurv, 3, 3);

  GaussCurvatureDEC *rhs = new GaussCurvatureDEC(sphere.getFeSpace(3));
  sphere.addVectorOperator(rhs,3);

  // ===== start adaption loop =====
  adapt->adapt();

  DOFVector<double> gcBonnet = *(sphere.getSolution(3));
  VtkVectorWriter::writeFile(gcBonnet, string("output/GaussBonnet.vtu"));

  DOFVector<double> mcMagY = halfMag(*(sphere.getSolution(0)), *(sphere.getSolution(1)), *(sphere.getSolution(2)), "laplaceMean", true);
  VtkVectorWriter::writeFile(mcMagY, "output/MeanMagY.vtu");

  MeshInfoCSVWriter mwriter("/dev/null/nonaynever.csv");
  mwriter.appendData(sphere.getFeSpace(),true);

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


