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


// Nonic Surface Pressed (double well (in z) strain of the unit sphere in x-direction (factor c on north pole and r*c on south pole)
// pressed to x-z-plane with factor 0<=b<1 (1:total flat pressed, 0:no pressing))
class PhiNP : public AbstractFunction<double, WorldVector<double> >
{
public:
  PhiNP(double c_,double r_,double b_) : AbstractFunction<double, WorldVector<double> >(1), c(c_), r(r_), b(b_) {}

  double operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    double dwell = (c*pow(z,2)*(r*(4 + 3*z)*pow(-1 + z,2) - (-4 + 3*z)*pow(1 + z,2)))/4.;
    double xcor = x - dwell;
    return xcor*xcor + y*y/((1.-b)*(1.-b)) + z*z - 1.0; 
  }

private:

  double c;
  double r;
  double b;
};

class GradPhiNP : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhiNP(double c_,double r_,double b_) : AbstractFunction<WorldVector<double>, WorldVector<double> >(1), c(c_), r(r_), b(b_) {}

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
    rval[1] = 2.0*y/((1.-b)*(1.-b));
    rval[2] = 2.0*(z - xcor*Ddwell);

    return rval;
  }

private:

  double c;
  double r;
  double b;
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
  TEST_EXIT(southRatio >= 0)("south radio factor must be positive!\n");

  double press = -1.0;
  Parameters::get("nonic->press", press);
  TEST_EXIT(press >= 0 && press < 1)("press factor must be positive and lower than 1!\n");

  double h;
  Parameters::get("meshCorrector->h", h);
  int nMax;
  Parameters::get("meshCorrector->nMax", nMax);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  string meshOut;
  Parameters::get("meshCorrector->outName", meshOut);
  MeshCorrector mc(sphere.getFeSpace(), meshOut);

// *********************************************
// all for nonicPressed
  double c;
  Parameters::get("edgeForces->c", c);
  double cNow = 0.8;
  Parameters::set("edgeForces->c", cNow);

  double hNow = h;

  double stretchOld = stretch;
  double pressOld = press;
  Parameters::get("nonic->old->c", stretchOld);
  Parameters::get("nonic->old->press", pressOld);

  int steps = 100;
  int iterationsPerStep = 101;
  Parameters::get("nonic->parameterSteps", steps);
  Parameters::get("nonic->iterationsPerStep", iterationsPerStep);


  double stretchStepSize = (stretchOld - stretch) / steps;
  double pressStepSize = (pressOld - press) / steps;

  double stretchNow = stretchOld - stretchStepSize;
  double pressNow = pressOld - pressStepSize;

  double signStretch = (stretchStepSize > 0 ) ? 1 : -1; 
  double signPress = (pressStepSize > 0 ) ? 1 : -1; 
  for (; signStretch*stretchNow > signStretch*stretch && signPress*pressNow > signPress*press; 
          stretchNow -= stretchStepSize, pressNow -= pressStepSize) {
    MSG("*** Precalculation at STRETCH = %f and PRESS = %f ***\n", stretchNow, pressNow);
    PhiProject proj(1, VOLUME_PROJECTION, new PhiNP(stretchNow, southRatio, pressNow), new GradPhiNP(stretchNow, southRatio, pressNow), 1.0e-6);
    mc.iterate(iterationsPerStep, hNow, meshOut);
    MSG("****************** DONE ***************************\n");
  }
  Parameters::set("edgeForces->c", c);
// *********************************************

  // ===== create projection =====
  //new PhiProject(1, VOLUME_PROJECTION, new PhiN(stretch, southRatio), new GradPhiN(stretch, southRatio), 1.0e-6);
  //new PhiProject(1, VOLUME_PROJECTION, new PhiO(stretch), new GradPhiO(stretch), 1.0e-6);
  new PhiProject(1, VOLUME_PROJECTION, new PhiNP(stretch, southRatio, press), new GradPhiNP(stretch, southRatio, press), 1.0e-9);
  //new TorusProject(1, VOLUME_PROJECTION, 2.0, 0.5);
  //WorldVector<double> ballCenter;
  //ballCenter.set(0.0);
  //new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);
  
  mc.iterate(nMax, h, meshOut);


  AMDiS::finalize();
}


