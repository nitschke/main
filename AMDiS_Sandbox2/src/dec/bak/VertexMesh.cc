#include "VertexMesh.h"

using namespace AMDiS;
using namespace dec;

VertexMesh::VertexMesh(const EdgeMesh *edgeMesh) : eMesh(edgeMesh) {
  FUNCNAME("VertexMesh::VertexMesh");
  MSG("Init VertexMesh ...\n");
  Timer t;

  int n = eMesh->getFeSpace()->getMesh()->getNumberOfAllDofs();

  vertices = vector<VertexElement>(n);

  vector<list<EdgeElement*> > edgesOnVertices(n);
  vector<EdgeElement>::iterator eIter = eMesh->getEdges()->begin();
  for(; eIter !=  eMesh->getEdges()->end(); ++eIter) {
    edgesOnVertices[eIter->dofEdge.first].push_back(&(*eIter));
    edgesOnVertices[eIter->dofEdge.second].push_back(&(*eIter));
  }

  MSG("done (needed %.5f seconds)\n", t.elapsed());
}
