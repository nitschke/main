#include "Dec.h"
#include "SphereProjection.h"
#include "phiProjection.h"
#include "EllipsoidProjection.h"
#include "torusProjection.h"
#include "ExtremeValueTracker.h"


using namespace std;
using namespace AMDiS;
using namespace dec;

// Rot(z)
class RotZ_Sphere : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotZ_Sphere() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x, const WorldVector<double>& vec) const 
  {
    double con2 = 1.0; //contra coords
    WorldVector<double> conBasis2; //contra basis vectors
    conBasis2[0] = x[1];
    conBasis2[1] = -x[0];
    conBasis2[2] = 0.0;
    return (con2 * conBasis2) * vec;
  }
};

class DZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = z
    return  q[2] - p[2];
  }
};

class DX : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DX() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = x
    return  q[0] - p[0];
  }
};


// Rot(z)
class RotXYZ_Sphere : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotXYZ_Sphere() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> hdfvec;
    hdfvec[0] = -0.25*x*(-1. + std::pow(x,2) - 3.*std::pow(y,2) + 5.*std::pow(z,2));
    hdfvec[1] = 0.25*y*(-1. - 3.*std::pow(x,2) + std::pow(y,2) + 5.*std::pow(z,2));
    hdfvec[2] = (std::pow(x,2) - 1.*std::pow(y,2))*z;
    return hdfvec * vec;
  }
};

class Xu_Torus : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  Xu_Torus() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldVector<double>operator()(const WorldVector<double>& coords) const 
  {
    WorldVector<double> coordss = coords;
    Projection::getProjection(1)->project(coordss);
    double x = coordss[0];
    double z = coordss[2];
    WorldVector<double> rvec;
    rvec[0] = -1.*z;
    rvec[1] = 0.;
    rvec[2] = x;
    return rvec;
  }
};

class XuHarm_Torus : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  XuHarm_Torus() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldVector<double>operator()(const WorldVector<double>& coords) const 
  {
    WorldVector<double> coordss = coords;
    Projection::getProjection(1)->project(coordss);
    double x = coordss[0];
    double z = coordss[2];
    WorldVector<double> rvec;
    rvec[0] = (-0.25*z)/(std::pow(x,2) + std::pow(z,2));
    rvec[1] = 0.;
    rvec[2] = (0.25*x)/(std::pow(x,2) + std::pow(z,2));
    return rvec;
  }
};

// Non-Harmonic part of Xu
class XuNonHarm_Torus : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  XuNonHarm_Torus() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldVector<double>operator()(const WorldVector<double>& coords) const 
  {
    WorldVector<double> coordss = coords;
    Projection::getProjection(1)->project(coordss);
    double x = coordss[0];
    double y = coordss[1];
    double z = coordss[2];
    WorldVector<double> rvec;
    rvec[0] = 0.125*z*(-4. + (47. + 4.*std::pow(y,2))/(std::pow(x,2) + std::pow(z,2)) - 16./std::sqrt(std::pow(x,2) + std::pow(z,2)));
    rvec[1] = 0.;
    rvec[2] = 0.125*x*(4. - (1.*(47. + 4.*std::pow(y,2)))/(std::pow(x,2) + std::pow(z,2)) + 16./std::sqrt(std::pow(x,2) + std::pow(z,2)));
    return rvec;
  }
};


// A_Xu*XuHarm + A_Xv*XvHarm
class LinHarm_Torus : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  LinHarm_Torus(double AXu_ = 0.5, double AXv_ = 0.5) : AbstractFunction<WorldVector<double>, WorldVector<double> >(), AXu(AXu_), AXv(AXv_) {}

  /// Implementation of AbstractFunction::operator().
  WorldVector<double>operator()(const WorldVector<double>& coords) const 
  {
    WorldVector<double> coordss = coords;
    Projection::getProjection(1)->project(coordss);
    double x = coordss[0];
    double y = coordss[1];
    double z = coordss[2];
    WorldVector<double> rvec;
    rvec[0] = AXu * (-0.25*z)/(std::pow(x,2) + std::pow(z,2)) + AXv * (-0.5*x*y)/(std::pow(x,2) + std::pow(z,2));
    rvec[1] = AXv * (0.5 - 1./std::sqrt(std::pow(x,2) + std::pow(z,2)));
    rvec[2] = AXu * (0.25*x)/(std::pow(x,2) + std::pow(z,2)) + AXv * (-0.5*y*z)/(std::pow(x,2) + std::pow(z,2));
    return rvec;
  }

private:
  double AXu;
  double AXv;
};

class DXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = xyz
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
  }
};

class DXYZZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    // f(X) = xyz
    return  q[0]*q[1]*q[2]*q[2] - p[0]*p[1]*p[2]*p[2];
  }
};

// d(Ax*X + Ay*Y + Az*Z)
class DLin : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DLin(double Ax_, double Ay_, double Az_) : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >(), Ax(Ax_), Ay(Ay_), Az(Az_) {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  Ax*(q[0]-p[0]) + Ay*(q[1]-p[1]) + Az*(q[2]-p[2]);
  }

private:
  double Ax;
  double Ay;
  double Az;
};


class GaussCurv_Sphere : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_Sphere() : AbstractFunction<double, EdgeElement >(){}

  double operator()(const EdgeElement& eel) const {
    return 1.0;
  }
};

class GaussCurv_Ellipsoid : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_Ellipsoid(double xscale, double yscale, double zscale) : AbstractFunction<double, EdgeElement >(), cx(xscale), cy(yscale), cz(zscale) {}

  double operator()(const EdgeElement& eel) const {
    WorldVector<double> coords = eel.infoLeft->getEdgeCenter(eel.dofEdge);
    Projection::getProjection(1)->project(coords);
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    double K = (4.*std::pow(cx,6)*std::pow(cy,6)*std::pow(cz,6))/std::pow(std::pow(cy,4)*std::pow(cz,4)*std::pow(x,2) + std::pow(cx,2)*std::pow(cy,2)*std::pow(cz,2)*(-1.*std::pow(cz,2)*(std::pow(x,2) + std::pow(y,2)) + std::pow(cy,2)*(std::pow(cz,2) - 1.*std::pow(z,2))) + std::pow(cx,4)*(std::pow(cz,4)*std::pow(y,2) + 2.*std::pow(cy,4)*std::pow(z,2) + std::pow(cy,2)*(std::pow(cz,4) - 1.*std::pow(cz,2)*std::pow(z,2))),2);    
    return K;
  }

  private:
  double cx; 
  double cy;
  double cz;
};

class GaussCurv_RBC : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_RBC(double a_,double c_) : AbstractFunction<double, EdgeElement >(), a(a_), c(c_) {
    FUNCNAME("GaussCurv_RBC::GaussCurv_RBC(double a_,double c_)")
    TEST_EXIT(a < c && c < a*std::sqrt(2))("It must hold: a < c < a*sqrt(2)");
  }

  double operator()(const EdgeElement& eel) const {
    WorldVector<double> coords = eel.infoLeft->getEdgeCenter(eel.dofEdge);
    Projection::getProjection(1)->project(coords);
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    double K = (-4.*std::pow(a,6) + std::pow(c,4)*std::sqrt(std::pow(c,4) + 4.*std::pow(a,2)*(std::pow(x,2) + std::pow(y,2))) + 8.*std::pow(a,4)*(-3.*std::pow(x,2) - 3.*std::pow(y,2) + std::sqrt(std::pow(c,4) + 4.*std::pow(a,2)*(std::pow(x,2) + std::pow(y,2)))) + std::pow(a,2)*(-5.*std::pow(c,4) + 6.*(std::pow(x,2) + std::pow(y,2))*std::sqrt(std::pow(c,4) + 4.*std::pow(a,2)*(std::pow(x,2) + std::pow(y,2)))))/(std::pow(c,4)*(std::pow(a,4) + std::pow(c,4) + std::pow(a,2)*(4.*std::pow(x,2) + 4.*std::pow(y,2) - 2.*std::sqrt(std::pow(c,4) + 4.*std::pow(a,2)*(std::pow(x,2) + std::pow(y,2))))));    
    return K;
  }

  private:
  double a;
  double c;
};

class GaussCurv_Torus : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_Torus(double r_ = 0.5,double R_ = 2.0) : AbstractFunction<double, EdgeElement >(), r(r_), R(R_) {}

  double operator()(const EdgeElement& eel) const {
    WorldVector<double> coords = eel.infoLeft->getEdgeCenter(eel.dofEdge);
    Projection::getProjection(1)->project(coords);
    double x = coords[0];
    double z = coords[2];
    return 4. - 8./std::sqrt(std::pow(x,2) + std::pow(z,2));
  }

  private:
  double r;
  double R;
};

class GaussCurv_Nonic095r : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_Nonic095r(double press, double stretch) : AbstractFunction<double, EdgeElement >(), c(stretch), B(press) {}

  double operator()(const EdgeElement& eel) const {
    WorldVector<double> coords = eel.infoLeft->getEdgeCenter(eel.dofEdge);
    Projection::getProjection(1)->project(coords);
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    double K = (-5.24288e11*std::pow(-1. + B,2)*(-3200. + 240.*c*x*(52. + 5.*z - 208.*std::pow(z,2) - 15.*std::pow(z,3) + 156.*std::pow(z,4) + 10.*std::pow(z,5)) + 3.*std::pow(c,2)*std::pow(z,2)*(-8112. - 1040.*z + 36479.*std::pow(z,2) + 3926.*std::pow(z,3) - 40470.*std::pow(z,4) - 4134.*std::pow(z,5) + 12073.*std::pow(z,6) + 1248.*std::pow(z,7) + 30.*std::pow(z,8))))/std::pow(4.096e7 + 1440.*std::pow(c,3)*x*std::pow(z,4)*std::pow(104. + 5.*z - 104.*std::pow(z,2) - 5.*std::pow(z,3),2)*(-156. - 5.*z + 78.*std::pow(z,2) + 3.*std::pow(z,3)) + 3.072e6*c*x*std::pow(z,2)*(-104. - 5.*z + 104.*std::pow(z,2) + 5.*std::pow(z,3)) + 9.*std::pow(c,4)*std::pow(z,6)*std::pow(16224. + 1300.*z - 24311.*std::pow(z,2) - 2002.*std::pow(z,3) + 8072.*std::pow(z,4) + 702.*std::pow(z,5) + 15.*std::pow(z,6),2) + 19200.*std::pow(c,2)*std::pow(z,2)*(-104. - 5.*z + 104.*std::pow(z,2) + 5.*std::pow(z,3))*(2.*std::pow(z,2)*(-156. - 5.*z + 78.*std::pow(z,2) + 3.*std::pow(z,3)) + 3.*std::pow(x,2)*(-104. - 5.*z + 104.*std::pow(z,2) + 5.*std::pow(z,3))) - 2.*B*(160.*c*x*std::pow(z,2)*(-2.9952e6 - 128000.*z - 9984.*(-250. + 1521.*std::pow(c,2))*std::pow(z,2) - 11520.*(-10. + 169.*std::pow(c,2))*std::pow(z,3) + 3.788226e7*std::pow(c,2)*std::pow(z,4) + 4.914747e6*std::pow(c,2)*std::pow(z,5) - 3.0161898e7*std::pow(c,2)*std::pow(z,6) - 3.988179e6*std::pow(c,2)*std::pow(z,7) + 7.419672e6*std::pow(c,2)*std::pow(z,8) + 1.019637e6*std::pow(c,2)*std::pow(z,9) + 45630.*std::pow(c,2)*std::pow(z,10) + 675.*std::pow(c,2)*std::pow(z,11)) + 6400.*std::pow(x,2)*(6400. + 9.*std::pow(c,2)*std::pow(z,2)*std::pow(104. + 5.*z - 104.*std::pow(z,2) - 5.*std::pow(z,3),2)) + std::pow(z,2)*(4.096e7 + 9.*std::pow(c,4)*std::pow(z,4)*std::pow(16224. + 1300.*z - 24311.*std::pow(z,2) - 2002.*std::pow(z,3) + 8072.*std::pow(z,4) + 702.*std::pow(z,5) + 15.*std::pow(z,6),2) + 6400.*std::pow(c,2)*std::pow(z,2)*(121680. + 9360.*z - 170177.*std::pow(z,2) - 13728.*std::pow(z,3) + 54486.*std::pow(z,4) + 4680.*std::pow(z,5) + 99.*std::pow(z,6)))) + std::pow(B,2)*(160.*c*x*std::pow(z,2)*(-2.9952e6 - 128000.*z - 9984.*(-250. + 1521.*std::pow(c,2))*std::pow(z,2) - 11520.*(-10. + 169.*std::pow(c,2))*std::pow(z,3) + 3.788226e7*std::pow(c,2)*std::pow(z,4) + 4.914747e6*std::pow(c,2)*std::pow(z,5) - 3.0161898e7*std::pow(c,2)*std::pow(z,6) - 3.988179e6*std::pow(c,2)*std::pow(z,7) + 7.419672e6*std::pow(c,2)*std::pow(z,8) + 1.019637e6*std::pow(c,2)*std::pow(z,9) + 45630.*std::pow(c,2)*std::pow(z,10) + 675.*std::pow(c,2)*std::pow(z,11)) + 6400.*std::pow(x,2)*(6400. + 9.*std::pow(c,2)*std::pow(z,2)*std::pow(104. + 5.*z - 104.*std::pow(z,2) - 5.*std::pow(z,3),2)) + std::pow(z,2)*(4.096e7 + 9.*std::pow(c,4)*std::pow(z,4)*std::pow(16224. + 1300.*z - 24311.*std::pow(z,2) - 2002.*std::pow(z,3) + 8072.*std::pow(z,4) + 702.*std::pow(z,5) + 15.*std::pow(z,6),2) + 6400.*std::pow(c,2)*std::pow(z,2)*(121680. + 9360.*z - 170177.*std::pow(z,2) - 13728.*std::pow(z,3) + 54486.*std::pow(z,4) + 4680.*std::pow(z,5) + 99.*std::pow(z,6)))),2);
    return K;
  }

  private:
  double c; //stretch
  double B; //press
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



class MyInstat : public DecProblemInstat {
public:
  MyInstat(DecProblemStat *probStat, DofEdgeVector initSolP, DofEdgeVector initSolD)
      : DecProblemInstat(probStat),
        solPrimal(initSolP),
        solDual(initSolD),
        tracker(probStat)
  {
    FUNCNAME("MyInstat::MyInstat(...)");
    
    string csvfn;
    Parameters::get(probStat->getName() + "->output->filename", csvfn);
    csvfn += "EKin.csv"; 

    csvout.open(csvfn.c_str(), ios::out);
    //csvout << "Time, Error" << endl;

    cout << setprecision(10);
    csvout << setprecision(10);

    double EKin = 0.5 * DofEdgeVectorPD::L2Norm2(solPrimal, solDual);
    csvout << 0.0 << "," << EKin << endl;
    cout << "### KinE: " << EKin << " ###" << endl;

    minusDivSolDual = -1.0 * solDual.divOnEdgeCenter();

    //DofEdgeVectorPD evecPD(solPrimal, solDual);
    //tracker.trackdownMinima(evecPD.getNormOnEdges(), 0.0);
  }


  void closeTimestep() {
    double time = t;
    DecProblemInstat::closeTimestep();
    solPrimal = statProb->getSolution(0);
    solDual =  statProb->getSolution(1);
    minusDivSolDual = - 1.0 * solDual.divOnEdgeCenter();

    double EKin = 0.5 * DofEdgeVectorPD::L2Norm2(solPrimal, solDual);
    csvout << time << "," << EKin << endl;
    cout << "### KinE: " << EKin << " ###" << endl;

    //DofEdgeVectorPD evecPD(solPrimal, solDual);
    //tracker.trackdownMinima(evecPD.getNormOnEdges(), time);
  }

  DofEdgeVector* getSolPrimal() {
    return &solPrimal;
  }

  DofEdgeVector* getSolDual() {
    return &solDual;
  }

  DofEdgeVector* getMinusDivSolDual() {
    return &minusDivSolDual;
  }

  ~MyInstat() {csvout.close();}

private:
  DofEdgeVector solPrimal;
  DofEdgeVector solDual;
  DofEdgeVector minusDivSolDual;

  OneMinValTrackerInPositiveZ tracker;

  ofstream csvout;
};

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  //SphereProject proj(42, VOLUME_PROJECTION);
  //new TorusProject(1, VOLUME_PROJECTION);
  double cx = 0.5;
  double cy = 0.5;
  double cz = 1.5;
  new EllipsoidProject(1, VOLUME_PROJECTION, cx, cy, cz);

  double nu = -1.0;
  Parameters::get("userParameter->kinematic_viscosity", nu);
  TEST_EXIT(nu >= 0.0)("kinematic_viscosity must be positive");

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);

  // Definition of alpha0 = [*dz, -dz]//
  DofEdgeVector alphaP(edgeMesh, "alphaPrimalInit");
  DofEdgeVector alphaD(edgeMesh, "alphaDualInit");
  //alphaP.interpol(new Xu_Torus());
  //alphaP.interpol(new XuNonHarm_Torus());
  //alphaP.interpol(new XuHarm_Torus());
  //alphaP.interpol(new XvHarm_Torus());
  //alphaP.interpol(new LinHarm_Torus(0.01, 0.99));
  //alphaD = alphaP.hodgeDual();
  //alphaP.interpolGL4(new RotXYZ_Sphere(), proj.getProjection(), proj.getJProjection());
  //alphaP.interpolGL4(new RotZ_Sphere(), proj.getProjection(), proj.getJProjection());
  //alphaD.set(new DXYZ());
  //alphaD.set(new DXYZZ());
  //alphaD.set(new DLin(0.1, 1.0, 0.1));
  alphaD.set(new DLin(0., 1.0, 0.1)); // for ellipsoid
  //alphaD.set(new DLin(0., 1., 1.)); // for RBC
  //alphaD *= -1.0;
  //alphaD.set(new DZ());
  //alphaD.set(new DX());
  alphaP = alphaD.hodgeDual();
  alphaD *= -1.0;


  MyInstat sphereInstat(&decSphere, alphaP, alphaD);

  // Gauss curvature on edge circumcenters
  DofEdgeVector K(edgeMesh, "K"); 
  //double press = 0.35;
  //double stretch = 1.0;
  //PhiProject proj(1, VOLUME_PROJECTION, new PhiNP(stretch, 0.95, press), new GradPhiNP(stretch, 0.95, press), 1.0e-6);
  //K.set(new GaussCurv_Nonic095r(press, stretch));
  K.set(new GaussCurv_Ellipsoid(cx, cy, cz));
  //K.set(new GaussCurv_Sphere());
  //K.set(new GaussCurv_Torus());

  //double a = 0.72;
  //double c = 0.75;
  //new PhiProject(1, VOLUME_PROJECTION, new PhiRBC(a, c), new GradPhiRBC(a, c), 1.0e-8);
  //K.set(new GaussCurv_RBC(a, c));
  K.writeFile("K.vtu");

// determine hodge dual
  EdgeOperator HodgeAlpha;
  HodgeAlpha.addTerm(new HodgeAtEdges());
  decSphere.addMatrixOperator(HodgeAlpha, 1, 0);

  EdgeOperator HAlpha;
  HAlpha.addTerm(new IdentityAtEdges(-1.0));
  decSphere.addMatrixOperator(HAlpha, 1, 1);

// time derivations
  EdgeOperator DtPrimal;
  DtPrimal.addTerm(new IdentityAtEdges());
  DtPrimal.setUhOld(alphaP, 0);
  decSphere.addMatrixOperator(DtPrimal, 0, 0, sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(DtPrimal, 0, sphereInstat.getInvTauPtr());

// diffusion
  EdgeOperator RotrotPrimal;
  RotrotPrimal.addTerm(new LaplaceBeltramiAtEdges(-nu));
  decSphere.addMatrixOperator(RotrotPrimal, 0, 0);

  EdgeOperator Gauss;
  Gauss.addTerm(new EdgeVecAtEdges(&K, NULL, -2.0*nu));
  decSphere.addMatrixOperator(Gauss,0,0);

// *alpha rot(alpha)   (convection)
  EdgeOperator HAlphaRotAlpha;
  HAlphaRotAlpha.addTerm(new RotAtEdgeCenterAndEdgeVecAtEdges(sphereInstat.getSolDual()));
  HAlphaRotAlpha.setUhOld(alphaP, 0);
  decSphere.addMatrixOperator(HAlphaRotAlpha, 0, 0);
  decSphere.addVectorOperator(HAlphaRotAlpha, 0);

  EdgeOperator DivHAlphaHAlpha;
  DivHAlphaHAlpha.addTerm(new EdgeVecAtEdges(sphereInstat.getMinusDivSolDual()));
  decSphere.addMatrixOperator(DivHAlphaHAlpha, 0, 1);



// dq  (convection + pressure)
  EdgeOperator ExDQ;
  ExDQ.addTerm(new ExteriorDerivativeAtEdges());
  decSphere.addMatrixOperator(ExDQ,0, 2);

// q = 0.5 <alpha,alpha> + p  (convection + pressure)
// linearize: <alpha,alpha> = 2 <alphaOld,alpha> - <alphaOld,alphaOld>
//    i.e.:  <alphaOld,alpha> + p - q =  0.5 <alphaOld,alphaOld>
  VertexOperator jAlphaOldAlpha;
  jAlphaOldAlpha.addTerm(new InterProdPartAtVertices(sphereInstat.getSolPrimal()));
  decSphere.addMatrixOperator(jAlphaOldAlpha, 2, 0);

  VertexOperator jHAlphaOldHAlpha;
  jHAlphaOldHAlpha.addTerm(new InterProdPartAtVertices(sphereInstat.getSolDual()));
  decSphere.addMatrixOperator(jHAlphaOldHAlpha, 2, 1);

  VertexOperator jAlphaOldAlphaOld;
  jAlphaOldAlphaOld.addTerm(new InterProdPartAtVertices(sphereInstat.getSolPrimal(), 0.5));
  jAlphaOldAlphaOld.setUhOldAtEdges(alphaP, 0);
  decSphere.addVectorOperator(jAlphaOldAlphaOld, 2);

  VertexOperator jHAlphaOldHAlphaOld;
  jHAlphaOldHAlphaOld.addTerm(new InterProdPartAtVertices(sphereInstat.getSolDual(), 0.5));
  jHAlphaOldHAlphaOld.setUhOldAtEdges(alphaD, 1);
  decSphere.addVectorOperator(jHAlphaOldHAlphaOld, 2);

  VertexOperator IdP;
  IdP.addTerm(new IdentityAtVertices());
  decSphere.addMatrixOperator(IdP, 2, 3);

  VertexOperator IdQ;
  IdQ.addTerm(new IdentityAtVertices(-1.0));
  decSphere.addMatrixOperator(IdP, 2, 2);

// conservation of mass
  VertexOperator DivAlpha;
  DivAlpha.addTerm(new DivAtVertices());
  decSphere.addMatrixOperator(DivAlpha, 3, 0);

// One Dof Condition to determine pressure
  decSphere.setValAtDof(3, 3, 0, 0, 0.0);

  //decSphere.assembleSystem();
  //cout << decSphere.getSysMat() << endl;
  //cout << decSphere.getRhs() << endl;

  sphereInstat.solve();

  AMDiS::finalize();
}
