#ifndef DEC_FWD_H
#define DEC_FWD_H

#include "AMDiS.h"

namespace AMDiS {
  
  namespace dec {
    class DecOperator;
    class DecOperatorTerm;
    class DecProblemStat;
    class DecProblemInstat;
    class DofEdgeVector;
    class EdgeMesh;
    class ElVolumesInfo2d;

    struct EdgeElement;

    typedef enum {
      UNDEFINEDSPACE = 0,
      VERTEXSPACE= 1,
      EDGESPACE = 2,
      FACESPACE = 3
    } SpaceType;

    typedef MTLTypes::MTLMatrix SparseMatrix;
    typedef MTLTypes::MTLVector DenseVector;

  }

}

#endif
