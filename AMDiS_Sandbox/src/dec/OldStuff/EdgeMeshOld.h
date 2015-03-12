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

  int getNumberOfEdges() const {return edgeRefs.size();} 
  
  DofEdge getEdge(int i) const {return edgeRefs[i];} 

  const vector<DofEdge>* getEdges() const {return &edgeRefs;}

private:
  const FiniteElemSpace *feSpace;
  vector<DofEdge> edgeRefs;

};

#endif
