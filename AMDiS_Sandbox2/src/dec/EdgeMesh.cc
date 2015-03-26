#include "EdgeMesh.h"

using namespace AMDiS;

EdgeMesh::EdgeMesh(const FiniteElemSpace *feSpace_): feSpace(feSpace_) {
  Flag fillFlags = Mesh::FILL_COORDS | Mesh::FILL_DET | Mesh::FILL_GRD_LAMBDA;

  Mesh *mesh = feSpace->getMesh();

  nEdges = mesh->getNumberOfEdges();

  edges = vector<EdgeElement>(nEdges);

  vector< pair<int,int> > macroIndex(nEdges); // <leftNeigh, rightNeigh>

  deque< MacroElement * > macroElements = mesh->getMacroElements();
  vector< ElVolumesInfo2d* > elVols(macroElements.size()); 
  int pos = 0;
  for (deque< MacroElement * >::iterator meIter = macroElements.begin(); 
              meIter != macroElements.end(); 
              ++meIter) {
    ElInfo2d *elInfo = new ElInfo2d(mesh);
    elInfo->setFillFlag(fillFlags);
    elInfo->fillMacroInfo(*meIter);
    elInfo->fillDetGrdLambda();
    elVols[(*meIter)->getIndex()] = new ElVolumesInfo2d(elInfo);
    //cout << (*meIter)->getIndex() <<endl;
    for (int i = 0; i < 3; i++) {
      //cout << "neigh " << i << " : " << (*meIter)->getNeighbour((i+2)%3)->getIndex() << endl;
      //need local oriented edge, ->getEdge() gives global oriented edges
      DofEdge dofEdge((*meIter)->getElement()->getDof(i,0), (*meIter)->getElement()->getDof((i+1)%3,0));
      if (dofEdge.first < dofEdge.second) {
        edges[pos].edgeDof = pos;
        edges[pos].dofEdge = dofEdge; //TODO: condition for boundaries

        macroIndex[pos].first = (*meIter)->getIndex();
        macroIndex[pos].second = (*meIter)->getNeighbour((i+2)%3)->getIndex();

        pos++;
      }
    }
  }

  for (int i = 0; i < edges.size(); ++i) {
    edges[i].infoLeft = elVols[macroIndex[i].first];
    edges[i].infoRight = elVols[macroIndex[i].second];
  }
}

DOFVector< list<EdgeElement> > EdgeMesh::getEdgeRings() const {
  DOFVector< list<EdgeElement> > edgeRings(feSpace, "EdgeRings");
  
  for (vector<EdgeElement>::const_iterator edgeIter = edges.begin();
       edgeIter != edges.end(); ++edgeIter) {
    edgeRings[edgeIter->dofEdge.first].push_back(*edgeIter);
    edgeRings[edgeIter->dofEdge.second].push_back(*edgeIter);
  }

  return edgeRings;
}


map<int, pair<ElVolumesInfo2d*, vector<EdgeElement> > > EdgeMesh::getFaceEdges() const {
  map<int, pair<ElVolumesInfo2d*, vector<EdgeElement> > > faceEdges;
  for (vector<EdgeElement>::const_iterator edgeIter = edges.begin();
       edgeIter != edges.end(); ++edgeIter) {
    int indexLeft = edgeIter->infoLeft->getElInfo()->getElement()->getIndex();
    faceEdges[indexLeft].first = edgeIter->infoLeft;   
    faceEdges[indexLeft].second.push_back(*edgeIter);   

    int indexRight = edgeIter->infoRight->getElInfo()->getElement()->getIndex();
    faceEdges[indexRight].first = edgeIter->infoRight;   
    faceEdges[indexRight].second.push_back(*edgeIter);   
  }
  return faceEdges;
}
