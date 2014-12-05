#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"
#include "torusProjection.h"
#include "meshCorrector.h"
#include "MeshHelper.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

// Heine 5.1(c)
//class Phi : public AbstractFunction<double, WorldVector<double> >
//{
//public:
//  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}
//
//  double operator()(const WorldVector<double>& x) const 
//  {
//    return 0.5 * (x[0]*x[0] + 4.0*x[1]*x[1] + (4.0/9.0)*x[2]*x[2] - 1.0);
//  }
//};
//
//class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
//{
//public:
//  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}
//
//  WorldVector<double> operator()(const WorldVector<double>& x) const 
//  {
//    WorldVector<double> rval(x);
//    rval[1] *= 4.0;
//    rval[2] *= 4.0/9.0;
//    return rval;
//  }
//};

//// Heine 5.1(b)
//class Phi : public AbstractFunction<double, WorldVector<double> >
//{
//public:
//  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}
//
//  double operator()(const WorldVector<double>& x) const 
//  {
//    double z2 = x[2] * x[2];
//    double xMz2 = x[0] - z2;
//    double yMz2 = x[1] - z2;
//    return 0.5 * (xMz2 * xMz2 + yMz2 * yMz2 + z2 - 1.0);
//  }
//};
//
//class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
//{
//public:
//  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}
//
//  WorldVector<double> operator()(const WorldVector<double>& x) const 
//  {
//    WorldVector<double> rval(x);
//    double z2 = x[2] * x[2];
//    rval[0] -= z2;
//    rval[1] -= z2;
//    rval[2] += 2.0 * x[2]* (2.0 * z2 - x[0] - x[1]);
//    return rval;
//  }
//};

//// Chmutov
//class Phi : public AbstractFunction<double, WorldVector<double> >
//{
//public:
//  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}
//
//  double operator()(const WorldVector<double>& x) const 
//  {
//    double x2 = x[0] * x[0];
//    double y2 = x[1] * x[1];
//    double z2 = x[2] * x[2];
//    double x4 = x2 * x2;
//    double y4 = y2 * y2;
//    double z4 = z2 * z2;
//    return x4 + y4 + z4 - x2 - y2 - z2 - 0.375;
//  }
//};
//
//class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
//{
//public:
//  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}
//
//  WorldVector<double> operator()(const WorldVector<double>& x) const 
//  {
//    WorldVector<double> rval(x);
//    rval *= -2.0;
//    for (int i = 0; i < 3; i++) rval[i] = rval[i] + 4.0 * x[i] * x[i] * x[i];
//    return rval;
//  }
//};



// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  //new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-6);
  //new TorusProject(1, VOLUME_PROJECTION, 2.0, 0.5);
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
  
  double h;
  Parameters::get("meshCorrector->h", h);
  int nMax;
  Parameters::get("meshCorrector->nMax", nMax);

  string meshOut;
  Parameters::get("meshCorrector->outName", meshOut);
  MeshCorrector mc(sphere.getFeSpace());
  mc.iterate(nMax, h, meshOut);


  // ===== start adaption loop =====
  adapt->adapt();


  //MacroWriter::writeMacro(new DataCollector<double>(sphere.getFeSpace()), meshOut.c_str());
  //DOFVector<double> *phi = new DOFVector<double>(sphere.getFeSpace(),"phi");
  //phi->interpol(new Phi());
  //VtkVectorWriter::writeFile(phi, meshOut + string("_phi.vtu"));

  AMDiS::finalize();
}


