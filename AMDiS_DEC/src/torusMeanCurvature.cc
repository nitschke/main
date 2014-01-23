#include "AMDiS.h"
#include "decOperator.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

class MeanCurv : public AbstractFunction<double, WorldVector<double> >
{
public:
  MeanCurv() : AbstractFunction<double, WorldVector<double> >(1) {}

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

  

  // ===== start adaption loop =====
  adapt->adapt();

  //cout << torus.getSystemMatrix(0,0)->getBaseMatrix() << endl;

  DOFVector<double> meanCurvDofVector(torus.getFeSpace(),"meanCurvature");
  meanCurvDofVector.interpol(new MeanCurv());
  VtkVectorWriter::writeFile(meanCurvDofVector, string("output/torusMeanCurvature.vtu"));

  torus.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


