#include "AMDiS.h"
#include <utility> //std::pair

using namespace AMDiS;
using namespace std;

class EdgeMesh {

public:
  EdgeMesh(const FiniteElemSpace *feSpace_);

  const FiniteElemSpace* getFeSpace() {return feSpace;}

  int getNumberOfEdges() {return edgeRefs.size();}
  
  DofEdge getEdge(int i) {return edgeRefs[i];} 

  vector<DofEdge>* getEdges() {return &edgeRefs;}

private:
  const FiniteElemSpace *feSpace;
  vector<DofEdge> edgeRefs;

};
