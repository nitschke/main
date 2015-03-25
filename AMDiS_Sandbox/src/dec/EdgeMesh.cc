#include "EdgeMesh.h"

using namespace AMDiS;

inline EdgeElement* getEdgeElementFromRow(vector<EdgeElement> edgeRow, DegreeOfFreedom vertex) {
  for(vector<EdgeElement>::iterator edgeRowIter = edgeRow.begin(); edgeRowIter != edgeRow.end(); ++edgeRowIter) {
    if ( edgeRowIter->nextVertex == vertex ) {
      return new EdgeElement(*edgeRowIter);
    }
  }
  cout << "irgendwas ist schief gelaufen" << endl;
}

EdgeMesh::EdgeMesh(const FiniteElemSpace *feSpace_): feSpace(feSpace_) {
  Mesh *mesh = feSpace->getMesh();
  nEdges = mesh->getNumberOfEdges();
  nVertices = mesh->getNumberOfVertices();

  vector<priority_queue<DegreeOfFreedom> > edgeQueues(nVertices);

  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL); el; el = stack.traverseNext(el)) {
    for (int i = 0; i < 3; i++) {
      edgeQueues[el->getElement()->getDof(i,0)].push(el->getElement()->getDof((i+1)%3,0));
      //cout << el->getElement()->getDof((i+1)%3,0) << endl;
    }
  }

  edges = vector<vector<EdgeElement> >(nVertices);

  //vector<priority_queue<DegreeOfFreedom> >::iterator qIter = edgeQueues.begin();
  //vector<vector<EdgeElement> >::iterator vIter = edges.begin();
  //for(; vIter != edges.end(); ++qIter, ++vIter){
  //  (*vIter) = vector<EdgeElement>(qIter->size());
  //  vector<EdgeElement>::iterator ventryIter = vIter->begin();
  //  for (; ventryIter != vIter->end(); ++ventryIter) {
  //    (*ventryIter).nextVertex = qIter->top();
  //    qIter->pop();
  //  }
  //}

  int index = 0;
  for(DegreeOfFreedom i = 0 ; i < nVertices; ++i) {
    edges[i] = vector<EdgeElement>(edgeQueues[i].size());
    for(vector<EdgeElement>::iterator edgeRowIter = edges[i].begin(); edgeRowIter != edges[i].end(); ++edgeRowIter) {
      edgeRowIter->nextVertex = edgeQueues[i].top();
      edgeQueues[i].pop();
      if (i < edgeRowIter->nextVertex) {
        edgeRowIter->index = index++;
      } else {
        edgeRowIter->companion = getEdgeElementFromRow(edges[edgeRowIter->nextVertex], i);
        edgeRowIter->index = edgeRowIter->companion->index;
        edgeRowIter->companion->companion = &(*edgeRowIter);
        //cout << edgeRowIter->index << " : ( " << i << " , " << edgeRowIter->nextVertex << " ) " << endl;
      }
    }
  }

 //for(DegreeOfFreedom i = 0 ; i < nVertices; ++i) {
 //   cout << i << endl;
 //   for(vector<EdgeElement>::iterator edgeRowIter = edges[i].begin(); edgeRowIter != edges[i].end(); ++edgeRowIter) {
 //     cout << edgeRowIter->index << " : ( " << i << " , " << edgeRowIter->nextVertex << " ) " << endl;
 //     //cout << "Companion: " << edgeRowIter->companion->index << " : ( " << edgeRowIter->companion->nextVertex << " , " << i << " ) " << endl;
 //   }
 // }


}

