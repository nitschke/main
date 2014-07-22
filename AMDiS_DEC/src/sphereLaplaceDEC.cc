#include "AMDiS.h"
#include "decOperator.h"
#include "MeshHelper.h"
#include "meshCorrector.h"
#include "DOFVHelper.h"

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
    //return -2.0 * x[0];
    return -6.0 * x[0] * x[2];
  }
};

class Sol : public AbstractFunction<double, WorldVector<double> >
{
public:
  Sol(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    //return x[0];
    return x[0] * x[2];
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
  sphere.setWriteAsmInfo(true);
  //sphere.setAssembleMatrixOnlyOnce(0, 0, false);


  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  // ===== create matrix operator =====
  //Operator matrixOperator(sphere.getFeSpace());
  //matrixOperator.addTerm(new Simple_SOT(-1.0));
  //sphere.addMatrixOperator(&matrixOperator, 0, 0);
  
  LBeltramiDEC decOperator(sphere.getFeSpace());
  sphere.addMatrixOperator(&decOperator, 0, 0);

  int degree = sphere.getFeSpace()->getBasisFcts()->getDegree();

  // ===== create rhs operator =====
  //Operator rhsOperator(sphere.getFeSpace());
  //rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  FunctionDEC rhsOperator(new F(degree), sphere.getFeSpace());

  sphere.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  DOFVector<double> solDOFV(sphere.getFeSpace(),"solDOFV");
  solDOFV.interpol(new Sol(0));
  VtkVectorWriter::writeFile(solDOFV, string("output/sol.vtu"));
  printError(*(sphere.getSolution(0)), solDOFV, "Error");

  //cout << sphere.getSystemMatrix(0,0)->getBaseMatrix() << endl;
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;
  //sphere.getSystemMatrix(0,0)->calculateNnz();
  //cout << "NNZ: " << sphere.getSystemMatrix(0,0)->getNnz() << endl;

  sphere.writeFiles(adaptInfo, true);

  //DOFVector<double> vorvol = getDualVols(sphere.getFeSpace());
  //cout << "Vol      : " << vorvol.sum() << endl;
  //cout << "Vol_Error: " << abs(vorvol.sum() - 4.0 * M_PI) << endl;
  //VtkVectorWriter::writeFile(vorvol, string("output/vorvol.vtu"));

  //DOFVector<double> radii = getVoronoiRadii(sphere.getFeSpace());
  //VtkVectorWriter::writeFile(radii, string("output/vorradii.vtu"));

  //DOFVector<int> nrOfCon = getConnections(sphere.getFeSpace());
  //VtkVectorWriter::writeFile(nrOfCon, string("output/connections.vtu"));

  //DOFVector<double> vol = get1RingVols(sphere.getFeSpace());
  //VtkVectorWriter::writeFile(vol, string("output/vol1Ring.vtu"));

  //DOFVector<WorldVector<double> > conForces = getConnectionForces(sphere.getFeSpace());
  //VtkVectorWriter::writeFile(conForces, string("output/conForces.vtu"));

  //MeshCorrector mc(sphere.getFeSpace());
  //int n = 1;
  //for (int i = 0; i < n; i++) {
  //  mc.oneIteration(0.1);
  //}
  //sphere.setFeSpace(mc.getFeSpace());
  //DOFVector<double> vol2 = get1RingVols(mc.getFeSpace());
  //VtkVectorWriter::writeFile(vol2, string("output/newVol1Ring.vtu"));


  AMDiS::finalize();
}


