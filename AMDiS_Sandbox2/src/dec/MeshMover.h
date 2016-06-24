#ifndef MESHMOVER_H
#define MESHMOVER_H

#include "Dec_fwd.h"

using namespace std;

namespace AMDiS { namespace dec {

//steal from MeshCorrector; why so much pointer?
class MeshMover {
public:
  
  MeshMover(EdgeMesh *eMesh, 
            BinaryAbstractFunction<WorldVector<double>, WorldVector<double>, double> *CoordsRefAndTimeToNewCoords);

  void move(double time);

private:
  FiniteElemSpace *feSpace;
  EdgeMesh *emesh;
  WorldVector<DOFVector<double> *  > coordsRef;
  WorldVector<DOFVector<double> *  > newCoords;
  BinaryAbstractFunction<WorldVector<double>, WorldVector<double>, double> *coordsFun;
};

}}

#endif
