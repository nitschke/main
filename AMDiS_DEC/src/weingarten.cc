#include "AMDiS.h"
#include "decOperator.h"
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
  

  int oh = 0;
  for (int i = 0; i < 3; i++) {
    // N
    sphere.addMatrixOperator(new SimpleDEC(sphere.getFeSpace(i+oh), sphere.getFeSpace(i+oh)), i+oh, i+oh);
    sphere.addVectorOperator(new DualPrimalNormalDEC(i, sphere.getFeSpace(i+oh)), i+oh);
    for (int j = 0; j < 3; j++) {
       int pos = matIndex(i,j) + oh + 3;
       // -II_ij
       SimpleDEC *II = new SimpleDEC(sphere.getFeSpace(pos), sphere.getFeSpace(pos));
       II->setFactor(-1.0);
       sphere.addMatrixOperator(II, pos, pos);
       // [d(N_j)]_i
       PrimalPrimalGradDEC *dN = new PrimalPrimalGradDEC(i, sphere.getFeSpace(pos), sphere.getFeSpace(j+oh));
       sphere.addMatrixOperator(dN, pos, j+oh);
    }
  }
  

  // ===== start adaption loop =====
  adapt->adapt();
  
  WorldMatrix<DOFVector<double> * > IIDV;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      IIDV[i][j] = sphere.getSolution(matIndex(i,j) + oh + 3);
    }
  }

  DOFVector<WorldVector<double> > eigDofVector = getEigenVals(IIDV);
  VtkVectorWriter::writeFile(eigDofVector, string("output/eigenVals.vtu"));

  DOFVector<double> k0 = getComp(0, eigDofVector, "k0");
  VtkVectorWriter::writeFile(k0, "output/k0.vtu");

  DOFVector<double> k1 = getComp(1, eigDofVector, "k1");
  VtkVectorWriter::writeFile(k1, "output/k1.vtu");


  DOFVector<double> gcWeingarten = prod01(eigDofVector);
  VtkVectorWriter::writeFile(gcWeingarten, "output/GaussWeingarten.vtu");
  DOFVector<double> gcWeingartenScaled = sqrtScale(gcWeingarten);
  VtkVectorWriter::writeFile(gcWeingartenScaled, "output/GaussWeingartenScaled.vtu");

  DOFVector<double> mcWeingarten = halfSum01(eigDofVector);
  VtkVectorWriter::writeFile(mcWeingarten, "output/MeanWeingarten.vtu");
  DOFVector<double> mcWeingartenScaled = sqrtScale(mcWeingarten);
  VtkVectorWriter::writeFile(mcWeingartenScaled, "output/MeanWeingartenScaled.vtu");

  MeshInfoCSVWriter mwriter("/dev/null/nonaynever.csv");
  mwriter.appendData(sphere.getFeSpace(),true);

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}



