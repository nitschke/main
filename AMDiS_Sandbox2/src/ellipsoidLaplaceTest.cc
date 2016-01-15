#include "Dec.h"
#include "io/VtkVectorWriter.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================
// LaplaceDeRham(SOL) + SOL
class RHS : public AbstractFunction<WorldVector<double>, WorldVector<double> > {
  public:
    RHS() : AbstractFunction<WorldVector<double>, WorldVector<double> >(){}

    WorldVector<double> operator()(const WorldVector<double>& coords) const {
      double x = coords[0];
      double y = coords[1];
      double z = coords[2];
      WorldVector<double> rval;

      rval[0] = (0.25*y*(59049.*(55937. - 1701.*pow(x,6) + 45.*pow(x,4)*(377. + 756.*pow(y,2)) + 12.*pow(y,2)*(8473. + 4524.*pow(y,2) + 1296.*pow(y,4)) - 9.*pow(x,2)*(8473. + 15080.*pow(y,2) + 9072.*pow(y,4))) - 2916.*(933403. + 100845.*pow(x,4) + 963864.*pow(y,2) + 322704.*pow(y,4) - 54.*pow(x,2)*(13387. + 14940.*pow(y,2)))*pow(z,2) + 50544.*(15149. - 5679.*pow(x,2) + 7572.*pow(y,2))*pow(z,4) - 7.881632e7*pow(z,6)))/pow(81.*(-5. + 3.*pow(x,2) - 12.*pow(y,2)) + 148.*pow(z,2),3);

      rval[1] = (0.0625*x*(59049.*(-11399. + 243.*pow(x,6) - 9.*pow(x,4)*(323. + 2268.*pow(y,2)) - 36.*pow(y,2)*(6373. + 6460.*pow(y,2) + 3024.*pow(y,4)) + 3.*pow(x,2)*(6373. + 38760.*pow(y,2) + 45360.*pow(y,4))) + 2916.*(209341. + 15795.*pow(x,4) + 1080.*pow(y,2)*(1699. + 1170.*pow(y,2)) - 90.*pow(x,2)*(1699. + 7020.*pow(y,2)))*pow(z,2) + 81648.*(-1987. + 601.*pow(x,2) - 7212.*pow(y,2))*pow(z,4) + 1.448288e7*pow(z,6)))/pow(81.*(-5. + 3.*pow(x,2) - 12.*pow(y,2)) + 148.*pow(z,2),3);

      rval[2] = (559872.*x*y*z*(-207. + 27.*pow(x,2) - 108.*pow(y,2) + 52.*pow(z,2)))/pow(81.*(-5. + 3.*pow(x,2) - 12.*pow(y,2)) + 148.*pow(z,2),3);

      return rval;
    }
};

class SOL : public AbstractFunction<WorldVector<double>, WorldVector<double> > {
  public:
    SOL() : AbstractFunction<WorldVector<double>, WorldVector<double> >(){}

    WorldVector<double> operator()(const WorldVector<double>& coords) const {
      double x = coords[0];
      double y = coords[1];
      double z = coords[2];
      WorldVector<double> rval;

      rval[0] = -2.*y;

      rval[1] = x/2.;

      rval[2] = 0.;

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
  //WorldVector<double> ballCenter;
  //ballCenter.set(0.0);
  //new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);
  
  DofEdgeVector rhs(edgeMesh, "RHS");
  rhs.interpol(new RHS());
  rhs.writeSharpFile("output/RHS.vtu", &sphere);

  DofEdgeVector sol(edgeMesh, "SOL");
  sol.interpol(new SOL());
  sol.writeSharpFile("output/SOL.vtu", &sphere);



// *** Edge Operators *** //
 // LaplaceDeRham
  double M1 = -1.0;
  EdgeOperator LaplaceOperator;
  LaplaceOperator.addTerm(new LaplaceBeltramiAtEdges());
  LaplaceOperator.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceOperator, 0, 0, &M1);

  EdgeOperator IOperator;
  IOperator.addTerm(new IdentityAtEdges());
  decSphere.addMatrixOperator(IOperator, 0, 0);


  EdgeOperator VOperator;
  VOperator.addTerm(new EdgeVecAtEdges(&rhs));
  decSphere.addVectorOperator(VOperator, 0);


  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();
  
  AMDiS::finalize();
}


