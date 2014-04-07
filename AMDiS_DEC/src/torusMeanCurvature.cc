#include "AMDiS.h"
#include "decOperator.h"
#include "MatrixHelper.h"
#include "MeshHelper.h"
#include "WorldVectorHelper.h"
#include "DOFVHelper.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

class MC : public AbstractFunction<double, WorldVector<double> >
{
public:
  MC() : AbstractFunction<double, WorldVector<double> >(1) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    double r = 0.5;
    double R = 2.0;
    double zeta =  (abs(x[1]) < r) ? (sqrt(r*r - x[1]*x[1])) : 0.0;
    if (x[0]*x[0] + x[2]*x[2] < R*R) zeta *= -1.0;
    double eta = zeta / (R + zeta);
    return 0.5 * (1.0 + eta) / r;
  }
};

class GC : public AbstractFunction<double, WorldVector<double> >
{
public:
  GC() : AbstractFunction<double, WorldVector<double> >(1) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    double r = 0.5;
    double R = 2.0;
    //double zeta =  (abs(x[1]) < r) ? (sqrt(r*r - x[1]*x[1])) : 0.0;
    double zeta =  (true) ? (sqrt(r*r - x[1]*x[1])) : 0.0;
    if (x[0]*x[0] + x[2]*x[2] < R*R) zeta *= -1.0;
    double eta = zeta / (R + zeta);
    return eta / (r * r);
  }
};

//Coords
class X_i : public AbstractFunction<double, WorldVector<double> >
{
public:
  X_i(int i_) : AbstractFunction<double, WorldVector<double> >(1), i(i_) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return x[i];
  }

protected:
  int i;
};

class TorusProject : public Projection
{
public:
  /// Constructor.
  TorusProject(int id, 
	ProjectionType type,
	double R_= 2.0,
  double r_= 0.5) 
    : Projection(id, type),
  R(R_),
  r(r_)
  {}

  /// Destructor.
  virtual ~TorusProject() {}

  /// Implementation of Projection::project();
  void project(WorldVector<double> &x) 
  {
    WorldVector<double> uTilde= x;
    uTilde[1]= 0.0;;
    WorldVector<double> u= uTilde*(R/norm(uTilde));
    WorldVector<double> xTilde= x - u;
    WorldVector<double> yTilde= xTilde*(r/norm(xTilde));
    x= u + yTilde;
  }

protected:
  double R;
  double r;
};


// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("torus main");

  AMDiS::init(argc, argv);

  new TorusProject(1, VOLUME_PROJECTION);

  // ===== create and init the scalar problem ===== 
  ProblemStat torus("torus");
  torus.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("torus->adapt", torus.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("torus->adapt",
					       &torus,
					       adaptInfo);
  
  for (int i = 0; i < 3; i++) {
    // ===== create matrix operator =====
    torus.addMatrixOperator(new SimpleDEC(torus.getFeSpace(i), torus.getFeSpace(i)), i, i);
    // ===== create rhs operator =====
    torus.addVectorOperator(new LBeltramiInteriorFunctionDEC(new X_i(i), torus.getFeSpace(i)), i);
  }
  
  // GaussCurv: geometric operator
  SimpleDEC *gaussCurv = new SimpleDEC(torus.getFeSpace(3),torus.getFeSpace(3) );
  torus.addMatrixOperator(gaussCurv, 3, 3);

  GaussCurvatureDEC *rhs = new GaussCurvatureDEC(torus.getFeSpace(3));
  torus.addVectorOperator(rhs,3);


  int oh = 4;
  for (int i = 0; i < 3; i++) {
    // N
    torus.addMatrixOperator(new SimpleDEC(torus.getFeSpace(i+oh), torus.getFeSpace(i+oh)), i+oh, i+oh);
    torus.addVectorOperator(new DualPrimalNormalDEC(i, torus.getFeSpace(i+oh)), i+oh);
    for (int j = 0; j < 3; j++) {
       int pos = matIndex(i,j) + oh + 3;
       // -II_ij
       SimpleDEC *II = new SimpleDEC(torus.getFeSpace(pos), torus.getFeSpace(pos));
       II->setFactor(-1.0);
       torus.addMatrixOperator(II, pos, pos);
       // [d(N_j)]_i
       PrimalPrimalGradDEC *dN = new PrimalPrimalGradDEC(i, torus.getFeSpace(pos), torus.getFeSpace(j+oh));
      torus.addMatrixOperator(dN, pos, j+oh);
    }
  }
  

  // ===== start adaption loop =====
  adapt->adapt();

  WorldMatrix<DOFVector<double> * > IIDV;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      IIDV[i][j] = torus.getSolution(matIndex(i,j) + oh + 3);
    }
  }

  DOFVector<WorldVector<double> > eigDofVector = getEigenVals(IIDV);
  VtkVectorWriter::writeFile(eigDofVector, string("output/eigenVals.vtu"));

  DOFVector<double> gcDOFV(torus.getFeSpace(),"GaussCurvExact");
  gcDOFV.interpol(new GC());
  VtkVectorWriter::writeFile(gcDOFV, string("output/gaussExact.vtu"));

  DOFVector<double> mcDOFV(torus.getFeSpace(),"MeanCurvExact");
  mcDOFV.interpol(new MC());
  VtkVectorWriter::writeFile(mcDOFV, string("output/meanExact.vtu"));

  DOFVector<double> gcBonnet = *(torus.getSolution(3));
  printError(gcBonnet, gcDOFV, "GaussBonnet");

  DOFVector<double> gcWeingarten = prod01(eigDofVector);
  printError(gcWeingarten, gcDOFV, "GaussWeingarten");

  DOFVector<double> mcMagY = halfMag(*(torus.getSolution(0)), *(torus.getSolution(1)), *(torus.getSolution(2)));
  printError(mcMagY, mcDOFV, "MeanMagY");

  DOFVector<double> mcWeingarten = halfSum01(eigDofVector);
  printError(mcWeingarten, mcDOFV, "MeanWeingarten");

  MeshInfoCSVWriter mwriter("/dev/null/nonaynever.csv");
  mwriter.appendData(torus.getFeSpace(),true);

  torus.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


