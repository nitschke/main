#include "Dec.h"
#include "SphereProjection.h"
#include "EllipsoidProjection.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// <df,[p,q]>
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

class DXYZ : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
  }
};

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

class RotXYZ_Ellipsoid : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotXYZ_Ellipsoid() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& coords, const WorldVector<double>& vec) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> hdfvec;
    hdfvec[0] = (9.*x*(-9. + 9.*std::pow(x,2) - 108.*std::pow(y,2) + 20.*std::pow(z,2))*std::sqrt((16.*std::pow(x,4)*std::pow(z,2) + 16.*std::pow(y,2)*(81. + 8.*(-9. + 2.*std::pow(y,2))*std::pow(z,2) + 16.*std::pow(z,4)) + std::pow(x,2)*(81. + 8.*(-9. + 16.*std::pow(y,2))*std::pow(z,2) + 16.*std::pow(z,4)))/std::pow(9. - 4.*std::pow(z,2),2)))/(-810. + 486.*std::pow(x,2) - 1944.*std::pow(y,2) + 296.*std::pow(z,2));
    hdfvec[1] = (9.*y*(9. + 27.*std::pow(x,2) - 36.*std::pow(y,2) - 20.*std::pow(z,2))*std::sqrt((16.*std::pow(x,4)*std::pow(z,2) + 16.*std::pow(y,2)*(81. + 8.*(-9. + 2.*std::pow(y,2))*std::pow(z,2) + 16.*std::pow(z,4)) + std::pow(x,2)*(81. + 8.*(-9. + 16.*std::pow(y,2))*std::pow(z,2) + 16.*std::pow(z,4)))/std::pow(9. - 4.*std::pow(z,2),2)))/(-810. + 486.*std::pow(x,2) - 1944.*std::pow(y,2) + 296.*std::pow(z,2));
    hdfvec[2] = (-162.*(std::pow(x,2) - 4.*std::pow(y,2))*z*std::sqrt((16.*std::pow(x,4)*std::pow(z,2) + 16.*std::pow(y,2)*(81. + 8.*(-9. + 2.*std::pow(y,2))*std::pow(z,2) + 16.*std::pow(z,4)) + std::pow(x,2)*(81. + 8.*(-9. + 16.*std::pow(y,2))*std::pow(z,2) + 16.*std::pow(z,4)))/std::pow(9. - 4.*std::pow(z,2),2)))/(-405. + 243.*std::pow(x,2) - 972.*std::pow(y,2) + 148.*std::pow(z,2));
    return  hdfvec * vec;
  }
};



int main(int argc, char* argv[])
{
  FUNCNAME("main");
  
  AMDiS::init(argc, argv);
  
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  //SphereProject sproj(42, VOLUME_PROJECTION);
  EllipsoidProject sproj(42, VOLUME_PROJECTION, 1.0, 0.5, 1.5);
  
  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);

// Definition of df //
  
  //BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dff = new DZ();
  //BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *rotff = new RotZ_Sphere();
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *dff = new DXYZ();
  BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *rotff = new RotXYZ_Ellipsoid();
  
  DofEdgeVector dfP(edgeMesh, "df");
  dfP.set(dff);
  dfP.writeSharpFile("output/dfP.vtu", &sphere);
  DofEdgeVector dfD(edgeMesh, "Rotf");
  dfD.interpolGL4(rotff, sproj.getProjection(), sproj.getJProjection());
  dfD.writeSharpFile("output/dfD.vtu", &sphere);
  dfD.writeFile("output/dfDForm.vtu");

  DofEdgeVector dfDh = dfP.hodgeDual();
  dfDh.writeSharpFile("output/dfDh.vtu", &sphere);
  dfDh.writeFile("output/dfDhForm.vtu");

  // Alpha = Alpha
  EdgeOperator Alpha;
  Alpha.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(Alpha, 0, 0);

  EdgeOperator AlphaRHS;
  AlphaRHS.addTerm(new EdgeVecAtEdges(&dfP));
  decSphere.addVectorOperator(AlphaRHS, 0);

  // *Alpha - *Alpha = 0;
  EdgeOperator HodgeAlpha;
  HodgeAlpha.addTerm(new HodgeAtEdges());
  decSphere.addMatrixOperator(HodgeAlpha, 1, 0);

  EdgeOperator HAlpha;
  HAlpha.addTerm(new IdentityAtEdges(-1.0));
  decSphere.addMatrixOperator(HAlpha, 1, 1);

  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();
  
  
  
  double hdia = edgeMesh->getMaxFaceDiameter();
  double hlen = edgeMesh->getMaxEdgeDiameter();
  cout << endl;
  cout << "max Diamter:  " << hdia << endl;
  cout << "max Length:   " << hlen << endl;

  double errL2Rel_hodge   = dfDh.errorL2Rel (dfD);
  double errMaxRel_hodge  = dfDh.errorMaxRel(dfD);
  cout << endl;
  cout << "hodge: RelError L2:  " <<  errL2Rel_hodge << endl;
  cout << "hodge: RelError Max: " << errMaxRel_hodge << endl;

  double errL2Rel_PSys   = decSphere.getEdgeSolution(0).errorL2Rel (dfP);
  double errMaxRel_PSys  = decSphere.getEdgeSolution(0).errorMaxRel(dfP);
  cout << endl;
  cout << "Sys Primal: RelError L2:  " <<  errL2Rel_PSys << endl;
  cout << "Sys Primal: RelError Max: " << errMaxRel_PSys << endl;

  double errL2Rel_DSys   = decSphere.getEdgeSolution(1).errorL2Rel (dfD);
  double errMaxRel_DSys  = decSphere.getEdgeSolution(1).errorMaxRel(dfD);
  cout << endl;
  cout << "Sys Dual: RelError L2:  " <<  errL2Rel_DSys << endl;
  cout << "Sys Dual: RelError Max: " << errMaxRel_DSys << endl;


  AMDiS::finalize();
}
