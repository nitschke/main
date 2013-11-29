#include "AMDiS.h"
#include "elVolumesInfo2d.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >
{
public:
  F(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return -2.0 * x[0];
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
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);



  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  //// === create adapt ===
  //AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
	//				       &sphere,
	//				       adaptInfo);
  
  //// ===== create matrix operator =====
  //Operator matrixOperator(sphere.getFeSpace());
  //matrixOperator.addTerm(new Simple_SOT);
  //sphere.addMatrixOperator(&matrixOperator, 0, 0);

  //// ===== create rhs operator =====
  //Operator rhsOperator(sphere.getFeSpace());

  //int degree = sphere.getFeSpace()->getBasisFcts()->getDegree();

  //rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  //sphere.addVectorOperator(&rhsOperator, 0);

  //// ===== start adaption loop =====
  //adapt->adapt();

  Mesh *mesh = sphere.getMesh();
  deque<MacroElement*>::iterator it;
  MSG("Nodes: %d\n", mesh->getNumberOfNodes());
  MSG("Leaves: %d\n", mesh->getNumberOfLeaves());
  MSG("All Elements: %d\n", mesh->getNumberOfElements());
  MSG("Macros: %d\n", mesh->getNumberOfMacros());

  TraverseStack stack;
  ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
  int elCounter = 0;
  while (elInfo) {
    ElVolumesInfo2d volInfo(elInfo);
    MSG("Element %d\n", ++elCounter);
    elInfo->getCoords().print();
    elInfo = stack.traverseNext(elInfo);
  }

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


