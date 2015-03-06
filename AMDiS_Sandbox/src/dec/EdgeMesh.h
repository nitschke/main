#ifndef EDGEMESH_H
#define EDGEMESH_H

#include "AMDiS.h"
#include <utility> //std::pair

using namespace AMDiS;
using namespace std;

class EdgeMesh {

public:
  EdgeMesh(const FiniteElemSpace *feSpace_);

  const FiniteElemSpace* getFeSpace() const {return feSpace;}


private:
  const FiniteElemSpace *feSpace;
  vector<vector<DegreeOfFreedom> > edges;

};

#endif
