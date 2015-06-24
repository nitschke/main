#include "EdgeMesh.h"

using namespace AMDiS;
using namespace dec;

EdgeMesh::EdgeMesh(const FiniteElemSpace *feSpace_): feSpace(feSpace_) {
  FUNCNAME("EdgeMesh::EdgeMesh");
  MSG("Init EdgeMesh ...\n");
  Timer t;
  
  Flag fillFlags = Mesh::FILL_COORDS | Mesh::FILL_DET | Mesh::FILL_GRD_LAMBDA;

  Mesh *mesh = feSpace->getMesh();
  //TODO: dofCompress in problemstat
  mesh->dofCompress();

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

  // ref volInfos
  for (int i = 0; i < edges.size(); ++i) {
    edges[i].infoLeft = elVols[macroIndex[i].first];
    edges[i].infoRight = elVols[macroIndex[i].second];
  }

  // get edge neighs TODO: improve
  // map macroIndex -> vector of the 3 edge facettes
  map<int, vector<EdgeElement*> > faceEdges;
  for (vector<EdgeElement>::iterator edgeIter = edges.begin();
       edgeIter != edges.end(); ++edgeIter) {
    faceEdges[edgeIter->infoLeft->getElInfo()->getElement()->getIndex()].push_back(&(*edgeIter));
    faceEdges[edgeIter->infoRight->getElInfo()->getElement()->getIndex()].push_back(&(*edgeIter));
  }
  // iterate over the macros
  map<int, vector<EdgeElement*> >::iterator feIter = faceEdges.begin();
  for (; feIter != faceEdges.end(); ++feIter) {
    Element *face = mesh->getMacroElement(feIter->first)->getElement();
    for (int i = 0; i < 3; ++i) {
      int ii = (i+1)%3;
      int iii = (i+2)%3;
      EdgeElement* e0 = (feIter->second)[i];
      EdgeElement* e1 = (feIter->second)[ii];
      EdgeElement* e2 = (feIter->second)[iii];
      // is e1 first or second?
      pair<EdgeElement*, EdgeElement*> edgeNeigh;
      if (e0->dofEdge.first == e1->dofEdge.first || e0->dofEdge.first == e1->dofEdge.second) {
        edgeNeigh.first = e1;
        edgeNeigh.second = e2;
      } else {
         edgeNeigh.first = e2;
         edgeNeigh.second = e1;
      }
      // is face left of e0?
      if (e0->infoLeft->getElInfo()->getElement() == face) {
        e0->edgesLeft = edgeNeigh;
      } else {
        e0->edgesRight = edgeNeigh;
      }
    }
  }

  MSG("done (%d Edges; needed %.5f seconds)\n", nEdges, t.elapsed());
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

double EdgeMesh::getVol() const {
  double vol = 0.0;
  for (vector<EdgeElement>::const_iterator edgeIter = edges.begin();
       edgeIter != edges.end(); ++edgeIter) {
    vol += edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);    
  }
  return vol;
}
