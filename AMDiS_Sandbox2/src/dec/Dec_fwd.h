#ifndef DEC_FWD_H
#define DEC_FWD_H

#include "AMDiS.h"

namespace AMDiS {
  
  namespace dec {
    class DofEdgeVector;
    class DofVertexVector;
    class AnimationWriter;
    class DecOperator;
    class DecOperatorTerm;
    class DecProblemStat;
    class DecProblemInstat;
    class EdgeMesh;
    class ElVolumesInfo2d;
    class SolverInterface;

    struct EdgeElement;

    typedef enum {
      UNDEFINEDSPACE = 0,
      VERTEXSPACE= 1,
      EDGESPACE = 2,
      FACESPACE = 3
    } SpaceType;

    typedef enum {
      FIRSTVERTEX = 1,
      SECONDVERTEX = 2,
      LEFTFACE = 3,
      RIGHTFACE = 4
    } EdgeRingIteratorType;

    //typedef  mtl::compressed2D<double> SparseMatrix;
    //typedef  mtl::dense_vector<double> DenseVector;
    typedef MTLTypes::MTLMatrix SparseMatrix;
    typedef MTLTypes::MTLVector DenseVector;

  }

}

#endif
