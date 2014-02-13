#include "AMDiS.h"
//#include <boost/numeric/mtl/mtl.hpp>

namespace AMDiS {

typedef mtl::matrix::dense2D<double> dMat;
typedef mtl::vector::dense_vector<double> dVec;

DOFVector<WorldVector<double> > getEigenVals(WorldMatrix< DOFVector<double> * > II) {
  DOFVector<WorldVector<double> > rval(II[0][0]->getFeSpace(), "EigenVals");
  DOFIterator<WorldVector<double> > rIter(&rval, USED_DOFS);
  rIter.reset();

  WorldMatrix< DOFIterator<double> * > IIIter;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      IIIter[i][j] = new DOFIterator<double>(II[i][j], USED_DOFS);
      IIIter[i][j]->reset();
    }
  }

  while(!rIter.end()) {
    dMat mat(3,3);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        mat[i][j] = *(*IIIter[i][j]);
      }
    }

    dVec eigs = mtl::matrix::qr_algo(mat, 200);

    for (int i = 0; i < 3; i++) {
      (*rIter)[i] = eigs[i];
    }
    
    rIter++;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        (*IIIter[i][j])++;
      }
    }
  }

  return rval;
}

}
