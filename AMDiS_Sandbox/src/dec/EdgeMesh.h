#ifndef EDGEMESH_H
#define EDGEMESH_H

#include "AMDiS.h"
#include <utility> //std::pair

using namespace AMDiS;
using namespace std;

struct EdgeElement {
  int index;  
  DegreeOfFreedom nextVertex;
  EdgeElement *companion;
  //edgeContainer *companion, *next(, pre);
  // ElInfo *faceLeft, *faceRight;
};

class EdgeMesh {

public:
  EdgeMesh(const FiniteElemSpace *feSpace_);

  const FiniteElemSpace* getFeSpace() const {return feSpace;}

  int getSize() const {return nEdges;}
  
  //all positiv oriented edges
  vector<DofEdge> getCompressedDofEdges() const {
    vector<DofEdge> dofEdges(nEdges);
    for (int i = 0 ; i <  nVertices; ++i) {
      for (int k = 0; edges[i][k].nextVertex > i && k < edges[i].size(); ++k) {
        dofEdges[edges[i][k].index] = DofEdge(i, edges[i][k].nextVertex);
      }
    }
    return dofEdges;
  }

  //TODO: no orientation yet
  vector<dofEdge> getEdgeRing const (DegreeOfFreedom i) {
    vector<DofEdge> dofEdges(edges[i].size());
    for (int k = 0; k < edges[i].size(); ++k) {
        dofEdges[edges[i][k].index] = DofEdge(i, edges[i][k].nextVertex);
    }
  }





private:
  const FiniteElemSpace *feSpace;
  vector<vector<EdgeElement> > edges;

  int nEdges; //compressed (in fact) number of edges (with companions -> 2*nEdges with no boundaries)
  int nVertices;
};

#endif
