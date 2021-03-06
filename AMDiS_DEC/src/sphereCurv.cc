#include "AMDiS.h"
#include "decOperator.h"
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


class GC : public AbstractFunction<double, WorldVector<double> >
{
public:
  GC() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& coord) const 
  {
    return 1.0;
  }
};

class MC : public AbstractFunction<double, WorldVector<double> >
{
public:
  MC() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& coord) const 
  {
    return 1.0;
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

  //MinusAngleDEC *ma = new MinusAngleDEC(sphere.getFeSpace(3));
  //sphere.addVectorOperator(ma,3);

  //SimplePrimalDEC *tp = new SimplePrimalDEC(sphere.getFeSpace(3));
  //tp->setFactor(2.0*M_PI);
  //sphere.addVectorOperator(tp,3);


  int oh = 4;
  //int oh = 7;
  for (int i = 0; i < 3; i++) {
    // N
    //sphere.addMatrixOperator(new SimpleDEC(sphere.getFeSpace(i+oh), sphere.getFeSpace(i+oh)), i+oh-3, i+oh-3);
    //sphere.addVectorOperator(new DualPrimalNormalDEC(i, sphere.getFeSpace(i+oh-3)), i+oh-3);
    for (int j = 0; j < 3; j++) {
       int pos = matIndex(i,j) + oh;
       // -II_ij
       SimpleDEC *II = new SimpleDEC(sphere.getFeSpace(pos), sphere.getFeSpace(pos));
       II->setFactor(-1.0);
       sphere.addMatrixOperator(II, pos, pos);
       // [d(N_j)]_i
       //PrimalPrimalGradDEC *dN = new PrimalPrimalGradDEC(i, sphere.getFeSpace(pos), sphere.getFeSpace(j+oh-3));
       //sphere.addMatrixOperator(dN, pos, j+oh-3);
       PrimalPrimalGradFunctionDEC *dN = new PrimalPrimalGradFunctionDEC(i, new X(j), sphere.getFeSpace(pos), sphere.getFeSpace(pos));
       dN->setFactor(-1.0);
       sphere.addVectorOperator(dN, pos);
    }
  }
  

  // ===== start adaption loop =====
  adapt->adapt();

  WorldMatrix<DOFVector<double> * > IIDV;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      IIDV[i][j] = sphere.getSolution(matIndex(i,j) + oh);
    }
  }

  DOFVector<WorldVector<double> > eigDofVector = getEigenVals(IIDV);
  VtkVectorWriter::writeFile(eigDofVector, string("output/eigenVals.vtu"));

  //DOFVector<double> avK = getAverage(getAverage(*(sphere.getSolution(3))));
  //VtkVectorWriter::writeFile(avK, string("output/GaussCurvAverage.vtu"));

  DOFVector<double> gcDOFV(sphere.getFeSpace(),"GaussCurvExact");
  gcDOFV.interpol(new GC());
  VtkVectorWriter::writeFile(gcDOFV, string("output/gaussExact.vtu"));

  DOFVector<double> mcDOFV(sphere.getFeSpace(),"MeanCurvExact");
  mcDOFV.interpol(new MC());
  VtkVectorWriter::writeFile(mcDOFV, string("output/meanExact.vtu"));

  DOFVector<double> gcBonnet = *(sphere.getSolution(3));
  printError(gcBonnet, gcDOFV, "GaussBonnet");
  VtkVectorWriter::writeFile(gcBonnet, string("output/GaussBonnet"));

  DOFVector<double> gcWeingarten = prod01(eigDofVector);
  printError(gcWeingarten, gcDOFV, "GaussWeingarten");
  VtkVectorWriter::writeFile(gcWeingarten, "output/GaussWeingarten.vtu");

  DOFVector<double> mcMagY = halfMag(*(sphere.getSolution(0)), *(sphere.getSolution(1)), *(sphere.getSolution(2)));
  printError(mcMagY, mcDOFV, "MeanMagY");
  VtkVectorWriter::writeFile(mcMagY, "output/MeanMagY.vtu");

  DOFVector<double> mcWeingarten = halfSum01(eigDofVector);
  printError(mcWeingarten, mcDOFV, "MeanWeingarten");
  VtkVectorWriter::writeFile(mcWeingarten, "output/MeanWeingarten.vtu");

  MeshInfoCSVWriter mwriter("/dev/null/nonaynever.csv");
  mwriter.appendData(sphere.getFeSpace(),true);

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


