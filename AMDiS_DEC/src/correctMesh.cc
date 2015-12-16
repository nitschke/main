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

// Octic Surface (double well (in z) strain of the unit sphere in x-direction)
class PhiO : public AbstractFunction<double, WorldVector<double> >
{
public:
  PhiO(double c_) : AbstractFunction<double, WorldVector<double> >(1), c(c_) {}

  double operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    double z2 = z*z;
    double xcor = x - c * z2 * (2.0 - z2);
    return xcor*xcor + y*y + z2 - 1.0; 
  }

private:

  double c;
};

class GradPhiO : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhiO(double c_) : AbstractFunction<WorldVector<double>, WorldVector<double> >(1), c(c_) {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    WorldVector<double> rval;
    double z2 = z*z;
    double xcor = x - c * z2 * (2.0 - z2);

    rval[0] = 2.0*xcor;
    rval[1] = 2.0*y;
    rval[2] = 2.0*(z + 4.0*c*z*(z2 - 1.0)*xcor);

    return rval;
  }

private:

  double c;
};

// Nonic Surface (double well (in z) strain of the unit sphere in x-direction (factor c on north pole and r*c on south pole)
class PhiN : public AbstractFunction<double, WorldVector<double> >
{
public:
  PhiN(double c_,double r_) : AbstractFunction<double, WorldVector<double> >(1), c(c_), r(r_) {}

  double operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    double dwell = (c*pow(z,2)*(r*(4 + 3*z)*pow(-1 + z,2) - (-4 + 3*z)*pow(1 + z,2)))/4.;
    double xcor = x - dwell;
    return xcor*xcor + y*y + z*z - 1.0; 
  }

private:

  double c;
  double r;
};

class GradPhiN : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhiN(double c_,double r_) : AbstractFunction<WorldVector<double>, WorldVector<double> >(1), c(c_), r(r_) {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    WorldVector<double> rval;
    double dwell = (c*pow(z,2)*(r*(4 + 3*z)*pow(-1 + z,2) - (-4 + 3*z)*pow(1 + z,2)))/4.;
    double xcor = x - dwell;
    double Ddwell = (c*z*(-8 - 15*z + r*(-8 + 15*z))*(-1 + pow(z,2)))/4.;

    rval[0] = 2.0*xcor;
    rval[1] = 2.0*y;
    rval[2] = 2.0*(z - xcor*Ddwell);

    return rval;
  }

private:

  double c;
  double r;
};



// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  double stretch = -1.0;
  Parameters::get("octic->c", stretch);
  TEST_EXIT(stretch >= 0)("stretch factor must be positive!\n");

  double southRatio = -1.0;
  Parameters::get("nonic->south ratio", southRatio);
  TEST_EXIT(southRatio >= 0)("stretch factor must be positive!\n");
  // ===== create projection =====
  //new PhiProject(1, VOLUME_PROJECTION, new PhiN(stretch, southRatio), new GradPhiN(stretch, southRatio), 1.0e-6);
  new PhiProject(1, VOLUME_PROJECTION, new PhiO(stretch), new GradPhiO(stretch), 1.0e-6);
  //new TorusProject(1, VOLUME_PROJECTION, 2.0, 0.5);
  //WorldVector<double> ballCenter;
  //ballCenter.set(0.0);
  //new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);
  

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


