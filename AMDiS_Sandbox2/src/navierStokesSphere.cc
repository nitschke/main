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


class MyInstat : public DecProblemInstat {
public:
  MyInstat(DecProblemStat *probStat, DofEdgeVector initSolP, DofEdgeVector initSolD)
      : DecProblemInstat(probStat),
        solPrimal(initSolP),
        solDual(initSolD)
  {
    FUNCNAME("MyInstat::MyInstat(...)");
    
    string csvfn;
    Parameters::get(probStat->getName() + "->output->filename", csvfn);
    csvfn += "RelKinErr.csv"; 

    csvout.open(csvfn.c_str(), ios::out);
    //csvout << "Time, Error" << endl;

    cout << setprecision(10);
    csvout << setprecision(10);

    double RelKinErr = std::abs(1.0 - (3.0/(8.0*M_PI)) * DofEdgeVectorPD::L2Norm2(solPrimal, solDual));
    csvout << 0.0 << "," << RelKinErr << endl;
    cout << "### RelKinErr: " << RelKinErr << " ###" << endl;
  }

  void closeTimestep() {
    double time = t;
    DecProblemInstat::closeTimestep();
    solPrimal = statProb->getSolution(0);
    solDual =  statProb->getSolution(1);

    double RelKinErr = std::abs(1.0 - (3.0/(8.0*M_PI)) * DofEdgeVectorPD::L2Norm2(solPrimal, solDual));
    csvout << time << "," << RelKinErr << endl;
    cout << "### RelKinErr: " << RelKinErr << " ###" << endl;
  }

  DofEdgeVector* getSolPrimal() {
    return &solPrimal;
  }

  DofEdgeVector* getSolDual() {
    return &solDual;
  }

  ~MyInstat() {csvout.close();}

private:
  DofEdgeVector solPrimal;
  DofEdgeVector solDual;

  ofstream csvout;
};

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  SphereProject proj(42, VOLUME_PROJECTION);

  double nu = -1.0;
  Parameters::get("userParameter->kinematic_viscosity", nu);
  TEST_EXIT(nu >= 0.0)("kinematic_viscosity must be positive");

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);

  // Definition of alpha0 = [*dz, -dz]//
  DofEdgeVector alphaP(edgeMesh, "alphaPrimalInit");
  DofEdgeVector alphaD(edgeMesh, "alphaDualInit");
  //alphaP.interpolGL4(new RotXYZ_Sphere(), proj.getProjection(), proj.getJProjection());
  alphaP.interpolGL4(new RotZ_Sphere(), proj.getProjection(), proj.getJProjection());
  //alphaD.set(new DXYZ());
  //alphaD *= -1.0;
  alphaD.set(new DZ());
  //alphaD.set(new DX());
  //alphaP = alphaD.hodgeDual();
  alphaD *= -1.0;

  MyInstat sphereInstat(&decSphere, alphaP, alphaD);

  // Gauss curvature on edge circumcenters
  DofEdgeVector K(edgeMesh, "K"); 
  K.set(new GaussCurv_Sphere());
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


  sphereInstat.solve();

  AMDiS::finalize();
}
