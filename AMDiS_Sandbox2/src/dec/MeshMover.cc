#include "MeshMover.h"
#include "EdgeMesh.h"

using namespace AMDiS;
using namespace dec;

MeshMover::MeshMover(EdgeMesh *eMesh, 
            BinaryAbstractFunction<WorldVector<double>, WorldVector<double>, double> *CoordsRefAndTimeToNewCoords) :
                emesh(eMesh), coordsFun(CoordsRefAndTimeToNewCoords) {

  feSpace = const_cast<FiniteElemSpace *>(emesh->getFeSpace());
  for (int i = 0; i < 3; i++) {
      coordsRef[i] = new DOFVector<double>(feSpace, "reference coord");
    }
    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int numBasFcts = basFcts->getNumber();
    std::vector<DegreeOfFreedom> localIndices(numBasFcts);
    DOFAdmin *admin = feSpace->getAdmin();
    DegreeOfFreedom dof;

    std::map<DegreeOfFreedom, bool> visited;
    TraverseStack stack;
    for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS); el; el = stack.traverseNext(el)) {
      basFcts->getLocalIndices(el->getElement(), admin, localIndices);
      for (int i = 0; i < 3; i++) {
        dof = localIndices[i];
        if (!visited[dof]) {
          for (int k = 0; k < 3; k++){
            (*(coordsRef[k]))[dof] = (el->getCoord(i))[k];
          }
          visited[dof] = true;
        }
      }
    }


    feSpace->getMesh()->setParametric(NULL);

  for (int i = 0; i < 3; i++) {
      (newCoords)[i] = new DOFVector<double>(feSpace, "new coord");
    }
}

void MeshMover::move(double time) {

  //WorldVector<DOFVector<double> * > *newCoords = new WorldVector<DOFVector<double> * >();
  //for (int i = 0; i < 3; i++) {
  //    (newCoords)[i] = new DOFVector<double>(feSpace, "new coord");
  //  }
  WorldVector<DOFIterator<double> * >  cIter;
  WorldVector<DOFIterator<double> * >  ncIter;
  for (int i = 0; i < 3; i++) {
    cIter[i] = new DOFIterator<double>(coordsRef[i], USED_DOFS);
    ncIter[i] = new DOFIterator<double>((newCoords)[i], USED_DOFS);
    cIter[i]->reset();
    ncIter[i]->reset();
  }

  while(!cIter[0]->end()) {
    WorldVector<double> coord, newCoord;
    for (int i = 0; i < 3; ++i) {
      coord[i] = *(*cIter[i]);
    }
    newCoord = (*coordsFun)(coord, time);
    for (int i = 0; i < 3; ++i) {
      *(*ncIter[i]) = newCoord[i];
    }

    for (int i = 0; i < 3; ++i) {
      ++(*cIter[i]);
      ++(*ncIter[i]);
    }
  }

  for (int i = 0; i < 3; i++) {
      delete cIter[i];
      delete ncIter[i];
  }

  Parametric *parametric = feSpace->getMesh()->getParametric(); 
  if (parametric) delete parametric;
  parametric = new ParametricFirstOrder(&newCoords);
  feSpace->getMesh()->setParametric(parametric);

  // how much is the fish?
  delete emesh;
  emesh = new EdgeMesh(feSpace); 
}
