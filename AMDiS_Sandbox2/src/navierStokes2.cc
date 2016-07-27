#include "Dec.h"
#include "SphereProjection.h"
#include "phiProjection.h"

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


class GaussCurv_Sphere : public AbstractFunction<double, EdgeElement > {
  public:
  GaussCurv_Sphere() : AbstractFunction<double, EdgeElement >(){}

  double operator()(const EdgeElement& eel) const {
    return 1.0;
  }
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


class MyInstat : public DecProblemInstat {
public:
  MyInstat(DecProblemStat *probStat, DofEdgeVector initSolP, DofEdgeVector initSolD)
      : DecProblemInstat(probStat),
        solPrimal(initSolP),
        solDual(initSolD) {}

  void closeTimestep() {
    DecProblemInstat::closeTimestep();
    solPrimal = statProb->getSolution(0);
    solDual =  statProb->getSolution(1);
  }

  DofEdgeVector* getSolPrimal() {
    return &solPrimal;
  }

  DofEdgeVector* getSolDual() {
    return &solDual;
  }

private:
  DofEdgeVector solPrimal;
  DofEdgeVector solDual;
};

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  //SphereProject sproj(42, VOLUME_PROJECTION);

  double nu = -1.0;
  Parameters::get("userParameter->kinematic_viscosity", nu);
  TEST_EXIT(nu >= 0.0)("kinematic_viscosity must be positive");

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);

  // Definition of alpha0 = [*dz, -dz]//
  DofEdgeVector alphaP(edgeMesh, "alphaPrimalInit");
  DofEdgeVector alphaD(edgeMesh, "alphaDualInit");
  //alphaP.interpolGL4(new RotXYZ_Sphere(), sproj.getProjection(), sproj.getJProjection());
  //alphaD.set(new DXYZ());
  //alphaD *= -1.0;
  //alphaD.set(new DZ());
  alphaD.set(new DX());
  alphaP = alphaD.hodgeDual();
  alphaD *= -1.0;

  MyInstat sphereInstat(&decSphere, alphaP, alphaD);

  // Gauss curvature on edge circumcenters
  double press = 0.35;
  double stretch = 1.0;
  PhiProject proj(1, VOLUME_PROJECTION, new PhiNP(stretch, 0.95, press), new GradPhiNP(stretch, 0.95, press), 1.0e-6);
  DofEdgeVector K(edgeMesh, "K"); 
  K.set(new GaussCurv_Nonic095r(press, stretch));
  K.writeFile("KNonic.vtu");

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
  decSphere.addMatrixOperator(HAlphaRotAlpha, 0, 0);

// dq  (convection + pressure)
  EdgeOperator ExDQ;
  ExDQ.addTerm(new ExteriorDerivativeAtEdges());
  decSphere.addMatrixOperator(ExDQ,0, 2);

// q = 0.5 <alpha,alpha> + p  (convection + pressure)
  VertexOperator jAlphaAlpha;
  jAlphaAlpha.addTerm(new InterProdPartAtVertices(sphereInstat.getSolPrimal(), 0.5));
  decSphere.addMatrixOperator(jAlphaAlpha, 2, 0);

  VertexOperator jHAlphaHAlpha;
  jHAlphaHAlpha.addTerm(new InterProdPartAtVertices(sphereInstat.getSolDual(), 0.5));
  decSphere.addMatrixOperator(jHAlphaHAlpha, 2, 1);

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
