#include "EdgeMesh.h"

using namespace AMDiS;

inline DofEdge orient(DofEdge edge);
inline bool contain(DofEdge edge, vector<DofEdge> edgeVec, int n); 

EdgeMesh::EdgeMesh(const FiniteElemSpace *feSpace_): feSpace(feSpace_) {
  Mesh *mesh = feSpace->getMesh();

  vector<priority_queue<DegreeOfFreedom> > edgeQueues(mesh->getNumberOfEdges());

  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL); el; el = stack.traverseNext(el)) {
    for (int i = 0; i < 3; i++) {
      edgeQueues[el->getElement()->getDof(i,0)].push(el->getElement()->getDof((i+1)%3,0));
    }
  }

  edges = vector<vector<DegreeOfFreedom> >(mesh->getNumberOfEdges());

  vector<priority_queue<DegreeOfFreedom> >::iterator qIter = edgeQueues.begin();
  vector<vector<DegreeOfFreedom> >::iterator vIter = edges.begin();
  for(; vIter != edges.end(); ++qIter, ++vIter){
    (*vIter) = vector<DegreeOfFreedom>(qIter->size());
    vector<DegreeOfFreedom>::iterator ventryIter = vIter->begin();
    for (; ventryIter != vIter->end(); ++ventryIter) {
      (*ventryIter) = qIter->top();
      qIter->pop();
    }
  }


}

