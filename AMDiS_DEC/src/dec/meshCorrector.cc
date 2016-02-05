#include "AMDiS.h"
#include "meshCorrector.h"
#include "MeshHelper.h"
#include "io/VtkVectorWriter.h"
#include "io/MacroWriter.h"


namespace AMDiS {


  MeshCorrector::MeshCorrector(const FiniteElemSpace *finiteElemSpace, std::string name) :
    feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace)),
    infowriter(string("meshStats" + name + ".csv"))
    {

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

    ii = 1;
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

  bool MeshCorrector::iterate(int n, double h, std::string name) {
    infowriter.appendData(feSpace);
    double tol1 = 1.0e-1;
    double tol2 = 1.0e-6;
    F = getConnectionForces(feSpace, true);
    if (ii == 1) {
      io::VtkVectorWriter::writeFile(F, string("output/ConForces" + name + "_" + boost::lexical_cast<std::string>(0) + ".vtu"));
    }

    double fMin = 1./0.;
    Parameters::get("corrector->fMin", fMin);

    int writeEveryIthStep = 10;
    Parameters::get("corrector->write every ith step", writeEveryIthStep);

    double fNew;
    double fOld = getMaxMagnitude(F);
    double hh = h;
    double hhOld = h;
    double fac = 1.2;
    double fac2 = 0.70;
    double k = 1.8;
    int nVerbose = 1;
    double h0 = 1.e-9;
    double n1 = 1000;
    int minusCounter = 0;
    int iMax = ii + n - 1;
    for (; ii < iMax || fOld > fMin; ii++) {
      //hh = (i > n1)? h : ((n1 - (double)i) * h0 + (double)i * h) / n1;
      //if (i%100 == 0 && i != 0) {
      //  fac = fac * fac * fac * fac * fac;
      //  fac2 = fac2 * fac2 * fac2 * fac2 * fac2;
      //}
      //oneHeunIteration(hh);
      oneIteration(hh);
      F = getConnectionForces(feSpace, true);
      fNew = getMaxMagnitude(F);
      //if (fNew < 1.0e-7) break;
      //hh *= 4.0 * fOld / (fNew + 3.0*fOld);
      double tmp = (fOld - fNew) / fOld;
      minusCounter = (tmp < 0) ? (minusCounter+1) : minusCounter;
      if (minusCounter > 100) {
        io::VtkVectorWriter::writeFile(F, string("output/ConForces" + name + "_EndAt_" + boost::lexical_cast<std::string>(ii) + ".vtu"));
        ERROR_EXIT("Mesh destroyed!");
        return false;
        hh *= 0.8;
        minusCounter = 0;
      }
      //if (tmp < tol1 && tmp > 0.0) {
      //  hh *= fac;
      //  fac2 = (1+k) - k * fac;
      //}
      //if (tmp < -tol2) {
      //    hh *= fac2;
      //    fac = 0.5 + 0.5 * fac;
      //}
      fOld = fNew;
      //hh *= (tmp < 0.01 && tmp >= 0) ? fac : fac2; 

      infowriter.appendData(feSpace);
      if (ii%nVerbose == 0) cout << ii << " : " << hh << " : " << fNew << " : " << tmp << " : " << fac << " : " << fac2 << endl;
      if (ii%writeEveryIthStep == 0) io::VtkVectorWriter::writeFile(F, string("output/ConForces" + name + "_" + boost::lexical_cast<std::string>(ii) + ".vtu"));
      if (ii%writeEveryIthStep == 0) {
        DataCollector<double> dc(feSpace);
        io::MacroWriter::writeMacro(&dc, string("output/meshOut" + name + "_" + boost::lexical_cast<std::string>(ii) + ".3d").c_str());
      }
    }

    return true;
  }


  WorldVector<DOFVector<double> * > MeshCorrector::coordsDeepCopy() {
    WorldVector<DOFVector<double> * > cccopy;
    for (int i = 0; i < 3; i++) {
      cccopy[i] = new DOFVector<double>(*(coords[i]));
    }
    return cccopy;
  }

}
