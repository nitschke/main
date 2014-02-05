#include "AMDiS.h"
#include "meshCorrector.h"
#include "MeshHelper.h"

namespace AMDiS {

  //// TODO: improve (calc coords only onetime, problem: coors must be set for all mels), 
  ////       projection from element without segfault
  //void MeshCorrector::oneIteration(double h) {
  //  DOFVector<WorldVector<double> > F = getConnectionForces(feSpace, true);

  //  TraverseStack stack;
  //  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS); el; el = stack.traverseNext(el)) {
  //    for (int i = 0; i < 3; i++) {
  //      DegreeOfFreedom dof = el->getMacroElement()->getElement()->getDof(i,0);
  //      //cout << el->getCoord(i) << endl;;
  //      WorldVector<double> newCoord = el->getMacroElement()->getCoord(i) + h * F[dof];
  //      //el->getMacroElement()->getProjection(i)->project(newCoord);
  //      Projection::getProjection(1)->project(newCoord);
  //      //cout << newCoord << endl;
  //      el->getMacroElement()->setCoord(i, newCoord);
  //      //cout << endl;
  //    }
  //  }
  //  //cout << "*************************" << endl;
  //}

  void MeshCorrector::oneIteration(double h) {
    DOFVector<WorldVector<double> > F = getConnectionForces(feSpace, true);

    WorldVector<DOFVector<double> * > newCoords;
    for (int i = 0; i < 3; i++) {
      newCoords[i] = new DOFVector<double>(feSpace, "paraCoordComponent");
    }

    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int numBasFcts = basFcts->getNumber();
    std::vector<DegreeOfFreedom> localIndices(numBasFcts);
    DOFAdmin *admin = feSpace->getAdmin();

    std::map<DegreeOfFreedom, bool> visited;
    TraverseStack stack;
    for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS); el; el = stack.traverseNext(el)) {
      basFcts->getLocalIndices(elInfo->getElement(), admin, localIndices);
      for (int i = 0; i < 3; i++) {
        dof = localIndices[i];
        if (!visited[dof]) {
          for (int k = 0; k < 3; k++){
            (*(newCoords[k]))[dof] = (el->getCoord(i))[k] + h * (F[dof])[k];
          }
          visited[dof] = true;
        }
      }
    }
    
  }

}
