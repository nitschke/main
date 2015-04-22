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


edgeRowValMapper LaplaceCoBeltramiAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  double f = fac * factor;

  // TODO: precalculatingn dual vols: ~O(6V)->O(V)
  // first dual volume
  double dvol1 = 0.0;
  for (EdgeElement::EdgeRingIterator eIter(&eel, FIRSTVERTEX); !eIter.isEnd(); ++eIter) {
    dvol1 += eIter.getFace()->getDualVertexVol_global(eel.dofEdge.first);
  }
  // second dual volume
  double dvol2 = 0.0;
  for (EdgeElement::EdgeRingIterator eIter(&eel, SECONDVERTEX); !eIter.isEnd(); ++eIter) {
    dvol2 += eIter.getFace()->getDualVertexVol_global(eel.dofEdge.second);
  }

  // first vertex iteration
  for (EdgeElement::EdgeRingIterator eIter(&eel, FIRSTVERTEX); !eIter.isEnd(); ++eIter) {
    // TODO: also precalculating: ~O(12E)->O(E)
    double metric = (  eIter->infoLeft->getDualEdgeLen(eIter->dofEdge)
                     + eIter->infoRight->getDualEdgeLen(eIter->dofEdge) )
                   / eIter->infoLeft->getEdgeLen(eIter->dofEdge);
    rowMapper[eIter->edgeDof] = f * ((eIter.pointOutward()) ? (-metric / dvol1) : (metric / dvol1));
  }
  //second vertex iteration
  EdgeElement::EdgeRingIterator eIter(&eel, SECONDVERTEX);
  // update on eel manually (first iterator object)
  rowMapper[eIter->edgeDof] -= f * ((  eIter->infoLeft->getDualEdgeLen(eIter->dofEdge)
                            + eIter->infoRight->getDualEdgeLen(eIter->dofEdge) )
                            / eIter->infoLeft->getEdgeLen(eIter->dofEdge)) / dvol2;
  for (++eIter; !eIter.isEnd(); ++eIter) {
    double metric = (  eIter->infoLeft->getDualEdgeLen(eIter->dofEdge)
                     + eIter->infoRight->getDualEdgeLen(eIter->dofEdge) )
                   / eIter->infoLeft->getEdgeLen(eIter->dofEdge);
    rowMapper[eIter->edgeDof] = f * ((eIter.pointOutward()) ? (metric / dvol2) : (-metric / dvol2));
  }

  return rowMapper;
}
