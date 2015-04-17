#include "EdgeOperatorTerm.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"
#include "ElVolumesInfo2d.h"

using namespace AMDiS;
using namespace dec;


edgeRowValMapper EdgeVecAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  DegreeOfFreedom edof = eel.edgeDof;
  rowMapper[edof] = fac * factor * (*evec)[edof];
  return rowMapper;
}

edgeRowValMapper IdentityAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  rowMapper[eel.edgeDof] = fac * factor;
  return rowMapper;
}


edgeRowValMapper LaplaceBeltramiAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  
  //duallength of eel (length of *eel)
  double dl = eel.infoLeft->getDualEdgeLen(eel.dofEdge); //restricted on left face
  dl += eel.infoRight->getDualEdgeLen(eel.dofEdge); //ristricted on right face

  // relation between legth and dual length
  double metric = eel.infoLeft->getEdgeLen(eel.dofEdge) / dl;

  //unsigned val on left face
  double valL = fac * factor * metric / eel.infoLeft->getVol();

  //unsigned val on right face
  double valR = fac * factor * metric / eel.infoRight->getVol();

  //left face iteration
  for (EdgeElement::EdgeRingIterator eIter(&eel, LEFTFACE); !eIter.isEnd(); ++eIter) {
    //respect orientation
    rowMapper[eIter->edgeDof] = (eIter.getType() == LEFTFACE) ? (-valL) : (valL);
  }

  //right face iteration
  EdgeElement::EdgeRingIterator eIter(&eel, RIGHTFACE);
  // update on eel manually (first iterator object)
  rowMapper[eIter->edgeDof] -= valR;
  for (++eIter; !eIter.isEnd(); ++eIter) {
    //respect orientation
    rowMapper[eIter->edgeDof] = (eIter.getType() == RIGHTFACE) ? (-valR) : (valR);
  }

  return rowMapper;
}


edgeRowValMapper LaplaceBeltramiAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;

  // TODO: precalculatingn dual vols: ~O(6V)->O(V)
  // first dual volume
  double dvol1 = 0.0;
  for (EdgeElement::EdgeRingIterator eIter(&eel, FIRSTVERTEX); !eIter.isEnd(); ++eIter) {
    dvol1 += eIter.getFace()->getDualVertexVol(eel->dofEdge.first);
  }
  // second dual volume
  double dvol2 = 0.0;
  for (EdgeElement::EdgeRingIterator eIter(&eel, SECONDVERTEX); !eIter.isEnd(); ++eIter) {
    dvol2 += eIter.getFace()->getDualVertexVol(eel->dofEdge.second);
  }

  // first vertex iteration
  for (EdgeElement::EdgeRingIterator eIter(&eel, FIRSTVERTEX); !eIter.isEnd(); ++eIter) {
    // TODO: also precalculating: ~O(12E)->O(E)
    double metric = (  eIter->InfoLeft->getDualEdgeLen(eIter->dofEdge)
                     + eIter->InfoRight->getDualEdgeLen(eIter->dofEdge) )
                   / eIter->InfoLeft->getEdgeLen(eIter->dofEdge);
    ...
  }

  return rowMapper;
}
