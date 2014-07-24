#include "AMDiS.h"
#include "decOperator.h"
#include "DOFVHelper.h"
#include "torusProjection.h"
#include "MeshHelper.h"

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
    //double zeta =  (abs(x[1]) < r) ? (sqrt(r*r - x[1]*x[1])) : 0.0;
    ////if (abs(x[1]) >= r) cout << x << endl;
    //if (x[0]*x[0] + x[2]*x[2] < R*R) zeta *= -1.0;
    ////if (abs(zeta) <= 0.01) cout << zeta << " on y = " << x[1] << endl; 
    //double eta = zeta / (R + zeta);
    //return - x[0] * eta * (eta + 1.0) / r / r;
    //return -x[1];
    //double x2Pz2 = x[0]*x[0] + x[2]*x[2];
    //return (x[0]/(r*r)) * ( (x[1]*x[1] - r*r)/x2Pz2 + R/sqrt(x2Pz2) - 1.0);
    double x2 = x[0]*x[0];
    double z2 = x[2]*x[2];
    //double x2Mz2 = x2 - z2;
    double x2Pz2 = x2 + z2;
    //return - 2.0 * x2Mz2 / (x2Pz2 * x2Pz2); 
    return 2.0 * z2 / (x2Pz2 * x2Pz2); 
  }
};

class G : public AbstractFunction<double, WorldVector<double> >
{
public:
  G(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    double x2 = x[0]*x[0];
    double z2 = x[2]*x[2];
    double x2Pz2 = x2 + z2;
    return 2.0 / x2Pz2; 
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
    double x2 = x[0]*x[0];
    double z2 = x[2]*x[2];
    double x2Pz2 = x2 + z2;
    return x2 / x2Pz2;
  }
};


// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("torus main");

  AMDiS::init(argc, argv);

  // ===== no projection, use finalize meshes =====
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
  
  // ===== create matrix operator =====
  //Operator laplaceOperator(torus.getFeSpace());
  //laplaceOperator.addTerm(new Simple_SOT(-1.0));
  LBeltramiDEC laplaceOperator(torus.getFeSpace());
  torus.addMatrixOperator(&laplaceOperator, 0, 0);

  int degree = torus.getFeSpace()->getBasisFcts()->getDegree();
  
  //Operator functionOperator(torus.getFeSpace());
  //functionOperator.addTerm(new CoordsAtQP_ZOT(new G(degree)));
  FunctionDEC functionOperator(new G(degree), torus.getFeSpace());
  torus.addMatrixOperator(&functionOperator, 0, 0);

  // ===== create rhs operator =====
  //Operator rhsOperator(torus.getFeSpace());
  //rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  FunctionDEC rhsOperator(new F(degree), torus.getFeSpace());
  torus.addVectorOperator(&rhsOperator, 0);


  

  // ===== start adaption loop =====
  adapt->adapt();

  //cout << torus.getSystemMatrix(0,0)->getBaseMatrix() << endl;

  DOFVector<double> rhsDofVector(torus.getFeSpace(),"RHS");
  rhsDofVector.interpol(new F(degree));
  VtkVectorWriter::writeFile(rhsDofVector, string("output/torusRHS.vtu"));

  DOFVector<double> solDOFV(torus.getFeSpace(),"solDOFV");
  solDOFV.interpol(new Sol(0));
  VtkVectorWriter::writeFile(solDOFV, string("output/sol.vtu"));
  printError(*(torus.getSolution(0)), solDOFV, "Error");

  MeshInfoCSVWriter mwriter("/dev/null/nonaynever.csv");
  mwriter.appendData(torus.getFeSpace(),true);  

  torus.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


