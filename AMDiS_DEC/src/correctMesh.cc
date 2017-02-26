#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"
#include "torusProjection.h"
#include "EllipsoidProjection.h"
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

//Red Blood Cell (RBC)
class PhiRBC : public AbstractFunction<double, WorldVector<double> >
{
public:
  PhiRBC(double a_,double c_) : AbstractFunction<double, WorldVector<double> >(1), a(a_), c(c_) {
    FUNCNAME("PhiRBC::PhiRBC(double a_,double c_)")
    TEST_EXIT(a < c && c < a*std::sqrt(2))("It must hold: a < c < a*sqrt(2)");
  }

  double operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    return -1.*std::pow(c,4) - 4.*std::pow(a,2)*(std::pow(x,2) + std::pow(y,2)) + std::pow(std::pow(a,2) + std::pow(x,2) + std::pow(y,2) + std::pow(z,2),2); 
  }

private:

  double a;
  double c;
};

class GradPhiRBC : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhiRBC(double a_,double c_) : AbstractFunction<WorldVector<double>, WorldVector<double> >(1), a(a_), c(c_) {
    FUNCNAME("GradPhiRBC::GradPhiRBC(double a_,double c_)")
    TEST_EXIT(a < c && c < a*std::sqrt(2))("It must hold: a < c < a*sqrt(2)");
  }

  WorldVector<double> operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    WorldVector<double> rval;
    
    rval[0] = 4.*x*(-1.*std::pow(a,2) + std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
    rval[1] = 4.*y*(-1.*std::pow(a,2) + std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
    rval[2] = 4.*z*(std::pow(a,2) + std::pow(x,2) + std::pow(y,2) + std::pow(z,2));

    return rval;
  }

private:

  double a;
  double c;
};



// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

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


  // ===== create projection =====
  //new PhiProject(1, VOLUME_PROJECTION, new PhiN(stretch, southRatio), new GradPhiN(stretch, southRatio), 1.0e-6);
  //new PhiProject(1, VOLUME_PROJECTION, new PhiO(stretch), new GradPhiO(stretch), 1.0e-6);
  //new PhiProject(1, VOLUME_PROJECTION, new PhiNP(stretch, southRatio, press), new GradPhiNP(stretch, southRatio, press), 1.0e-9);
  //new TorusProject(1, VOLUME_PROJECTION, 2.0, 0.5);
  //new EllipsoidProject(1, VOLUME_PROJECTION, 1.0, 1.0, 1.25);
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  //double a = 0.72;
  //double c = 0.75;
  //new PhiProject(1, VOLUME_PROJECTION, new PhiRBC(a, c), new GradPhiRBC(a, c), 1.0e-8);
  
  mc.iterate(nMax, h, meshOut);


  AMDiS::finalize();
}


