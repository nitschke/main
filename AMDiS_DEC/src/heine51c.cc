#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"
#include "meshCorrector.h"
#include "MeshHelper.h"

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

class Phi : public AbstractFunction<double, WorldVector<double> >
{
public:
  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& x) const 
  {
    return 0.5 * (x[0]*x[0] + 4.0*x[1]*x[1] + (4.0/9.0)*x[2]*x[2] - 1.0);
  }
};

class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}

  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> rval(x);
    rval[1] *= 4.0;
    rval[2] *= 4.0/9.0;
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
  new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-6);

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

  DOFVector<WorldVector<double> > forces = getConnectionForces(sphere.getFeSpace());
  VtkVectorWriter::writeFile(forces, string("output/ConForces.vtu"));

  DOFVector<WorldVector<double> > normals = getNormals(sphere.getFeSpace());
  VtkVectorWriter::writeFile(normals, string("output/Normals.vtu"));

  DOFVector<double> radii;

  double h;
  Parameters::get("meshCorrector->h", h);
  int nMax;
  Parameters::get("meshCorrector->nMax", nMax);

  MeshCorrector mc(sphere.getFeSpace());
  mc.iterate(nMax, h);
  //for (int i = 0; i < nMax; i++) {
  //  mc.iterate(1, h);
  //  
  //  forces = getConnectionForces(sphere.getFeSpace(), true);
  //  VtkVectorWriter::writeFile(forces, string("output/ConForces_" + boost::lexical_cast<std::string>(i) + ".vtu"));

  //  radii = getVoronoiRadii(sphere.getFeSpace());
  //  radii = radii.average();
  //  VtkVectorWriter::writeFile(radii, string("output/Radii_" + boost::lexical_cast<std::string>(i) + ".vtu"));
  //}

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


