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
  //

  MeshCorrector::MeshCorrector(const FiniteElemSpace *finiteElemSpace) :
    feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace)) {

    for (int i = 0; i < 3; i++) {
      coords[i] = new DOFVector<double>(feSpace, "paraCoordComponent");
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
            (*(coords[k]))[dof] = (el->getCoord(i))[k];
          }
          visited[dof] = true;
        }
      }
    }


    feSpace->getMesh()->setParametric(NULL);
  }

  void MeshCorrector::oneIteration(double h) {
    DOFIterator<WorldVector<double> > fIter(&F, USED_DOFS);
    fIter.reset();

    WorldVector<DOFIterator<double> * >  cIter;
    for (int i = 0; i < 3; i++) {
      cIter[i] = new DOFIterator<double>(coords[i], USED_DOFS);
      cIter[i]->reset();
    }
    
    while(!fIter.end()){
      WorldVector<double> newCoord;

      for (int i = 0; i < 3; i++) newCoord[i] = *(*(cIter[i])) + h * (*fIter)[i];
      Projection::getProjection(1)->project(newCoord);
      for (int i = 0; i < 3; (*(cIter[i]))++, i++) *(*(cIter[i])) = newCoord[i];

      fIter++;
    }

    for (int i = 0; i < 3; i++) {
      delete cIter[i];
    }

    Parametric *parametric = feSpace->getMesh()->getParametric(); 
    if (parametric) delete parametric;
    parametric = new ParametricFirstOrder(&coords);
    feSpace->getMesh()->setParametric(parametric);
  }

  void MeshCorrector::oneHeunIteration(double h) {
    WorldVector<DOFVector<double> * > oldCoords = coordsDeepCopy();
    oneIteration(h);
    F = getConnectionForces(feSpace, true);
    oneIteration(h);

    WorldVector<DOFIterator<double> * >  cIter;
    WorldVector<DOFIterator<double> * >  cOldIter;
    for (int i = 0; i < 3; i++) {
      cIter[i] = new DOFIterator<double>(coords[i], USED_DOFS);
      cIter[i]->reset();
      cOldIter[i] = new DOFIterator<double>(oldCoords[i], USED_DOFS);
      cOldIter[i]->reset();
    }

    while(!cIter[0]->end()){
      WorldVector<double> newCoord;
      for (int i = 0; i < 3;(*(cOldIter[i]))++, i++) newCoord[i] = 0.5 * ( *(*(cIter[i])) + *(*(cOldIter[i])));
      Projection::getProjection(1)->project(newCoord);
      for (int i = 0; i < 3; (*(cIter[i]))++, i++) *(*(cIter[i])) = newCoord[i];
    }

    for (int i = 0; i < 3; i++) {
      delete cIter[i];
      delete cOldIter[i];
    }

    Parametric *parametric = feSpace->getMesh()->getParametric(); 
    if (parametric) delete parametric;
    parametric = new ParametricFirstOrder(&coords);
    feSpace->getMesh()->setParametric(parametric);
  }

  void MeshCorrector::iterate(int n, double h) {
    double tol1 = 1.0e-1;
    double tol2 = 1.0e-6;
    F = getConnectionForces(feSpace, true);
    double fNew;
    double fOld = getMaxMagnitude(F);
    double hh = h;
    double hhOld = h;
    double fac = 1.5;
    double fac2 = 0.75;
    double k = 1.8;
    for (int i = 0; i < n; i++) {
      //if (i%100 == 0 && i != 0) {
      //  fac = fac * fac * fac * fac * fac;
      //  fac2 = fac2 * fac2 * fac2 * fac2 * fac2;
      //}
      oneHeunIteration(hh);
      F = getConnectionForces(feSpace, true);
      fNew = getMaxMagnitude(F);
      //hh *= 4.0 * fOld / (fNew + 3.0*fOld);
      double tmp = (fOld - fNew) / fOld;
      if (tmp < tol1 && tmp > 0.0) {
        hh *= fac;
        fac2 = (1+k) - k * fac;
      }
      if (tmp < -tol2) {
          hh *= fac2;
          fac = 0.5 + 0.5 * fac;
      }
      fOld = fNew;

      cout << i << " : " << hh << " : " << fNew << " : " << tmp << " : " << fac << " : " << fac2 << endl;
      VtkVectorWriter::writeFile(F, string("output/ConForces_" + boost::lexical_cast<std::string>(i) + ".vtu"));
    }
  }


  WorldVector<DOFVector<double> * > MeshCorrector::coordsDeepCopy() {
    WorldVector<DOFVector<double> * > cccopy;
    for (int i = 0; i < 3; i++) {
      cccopy[i] = new DOFVector<double>(*(coords[i]));
    }
    return cccopy;
  }

}
