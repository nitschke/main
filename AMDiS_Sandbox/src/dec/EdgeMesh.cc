#include "EdgeMesh.h"

using namespace AMDiS;

inline DofEdge orient(DofEdge edge);
inline bool contain(DofEdge edge, vector<DofEdge> edgeVec, int n); 

EdgeMesh::EdgeMesh(const FiniteElemSpace *feSpace_): feSpace(feSpace_) {
  Mesh *mesh = feSpace->getMesh();

  edgeRefs = vector<DofEdge>(mesh->getNumberOfEdges());

  int pos = 0;
  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL); el; el = stack.traverseNext(el)) {
    for (int i = 0; i < 3; i++) {
      DofEdge edge = orient(el->getElement()->getEdge(i));
      if (!contain(edge, edgeRefs, pos-1)) edgeRefs[pos++] = edge;
    }
  }

  //for (vector<DofEdge>::iterator edgeIter = edgeRefs.begin(); edgeIter != edgeRefs.end(); ++edgeIter)
  //  cout << (*edgeIter).first << " " << (*edgeIter).second << endl;


}

inline DofEdge orient(DofEdge edge) {
  return (edge.first > edge.second) ? DofEdge(edge.second,edge.first) : DofEdge(edge);
}

inline bool contain(DofEdge edge, vector<DofEdge> edgeVec, int n) {
  for (int i = 0; i <= n; i++) 
    if (edgeVec[i] == edge) return true;
  return false;
}
