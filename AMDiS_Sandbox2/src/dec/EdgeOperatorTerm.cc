#include "EdgeOperatorTerm.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"
#include "ElVolumesInfo2d.h"

using namespace AMDiS;
using namespace dec;


edgeRowValMapper EdgeVecAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  double f = fac * factor;
  double val = (*evec)[eel.edgeDof];
  rowMapper[eel.edgeDof] = f * ((func) ? (*func)(val) : val);
  return rowMapper;
}

edgeRowValMapper EdgeFunAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  rowMapper[eel.edgeDof] = fac * factor * (*func)(eel);
  return rowMapper;
}


edgeRowValMapper EdgeVec2AndEdgeAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  rowMapper[eel.edgeDof] = fac * factor * (*func)((*evec1)[eel], (*evec2)[eel], eel);
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

edgeRowValMapper NormSquaredEdgeVecAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  double f = fac * factor;

  ElVolumesInfo2d *infoL = eel.infoLeft;
  ElVolumesInfo2d *infoR = eel.infoRight;

  double detL = infoL->getElInfo()->getDet();
  double detR = infoR->getElInfo()->getDet();
  double det2InvL = 1.0 / (detL * detL);
  double det2InvR = 1.0 / (detR * detR);

  double edgeVal = (*evec)[eel];

  WorldVector<double> edge = infoL->getEdge(eel.dofEdge);

  EdgeElement *eelAdj;
  WorldVector<double> edgeAdj;
  WorldVector<double> evalVec;

  // left first
  eelAdj = eel.edgesLeft.first;
  edgeAdj = infoL->getEdge(eelAdj->dofEdge);
  evalVec = edgeVal * edgeAdj - (*evec)[eelAdj] * edge;
  double scaleL1 = infoL->getDualVertexVol_global(eel.dofEdge.first);
  double norm2L1 = det2InvL * (evalVec * evalVec);

  // left second
  eelAdj = eel.edgesLeft.second;
  edgeAdj = infoL->getEdge(eelAdj->dofEdge);
  evalVec = edgeVal * edgeAdj - (*evec)[eelAdj] * edge;
  double scaleL2 = infoL->getDualVertexVol_global(eel.dofEdge.second);
  double norm2L2 = det2InvL * (evalVec * evalVec);
 

  // right first
  eelAdj = eel.edgesRight.first;
  edgeAdj = infoR->getEdge(eelAdj->dofEdge);
  evalVec = edgeVal * edgeAdj - (*evec)[eelAdj] * edge;
  double scaleR1 = infoR->getDualVertexVol_global(eel.dofEdge.first);
  double norm2R1 = det2InvR * (evalVec * evalVec);

  // right second
  eelAdj = eel.edgesRight.second;
  edgeAdj = infoR->getEdge(eelAdj->dofEdge);
  evalVec = edgeVal * edgeAdj - (*evec)[eelAdj] * edge;
  double scaleR2 = infoR->getDualVertexVol_global(eel.dofEdge.second);
  double norm2R2 = det2InvR * (evalVec * evalVec);

  //averaging
  double rval = (scaleL1*norm2L1 + scaleL2*norm2L2 + scaleR1*norm2R1 + scaleR2*norm2R2) 
                              / (scaleL1 + scaleL2 + scaleR1 + scaleR2);
  rowMapper[eel.edgeDof] = f * ((func) ? (*func)(rval) : rval);

  return rowMapper;
}

//<dh,e> = h(v_2) - h(v_1)
edgeRowValMapper ExteriorDerivativeAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  double f = fac * factor;
  rowMapper[eel.dofEdge.second] = f;
  rowMapper[eel.dofEdge.first] = -f;

  return rowMapper;
}


/*            w2 
 *           /^\
 *          / | \
 *         / h|  \
 *       v1------>v2
 *         \ e|  /
 *          \ | /
 *           \|/
 *            w1
 */
// <*df,e> = (1 / sqrt(|g|))(<df,e>(e.h) - <df,h>|e|^2)
// |g| = |e|^2 |h|^2 - (e.h)^2
// e = [v1,v2] ; h = [w1 , w2]
edgeRowValMapper RotAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  double f = fac * factor;
  
  DegreeOfFreedom dofv1 = eel.dofEdge.first;
  DegreeOfFreedom dofv2 = eel.dofEdge.second;
  DegreeOfFreedom dofw2 = eel.infoLeft->getOppVertexDof (eel.dofEdge);
  DegreeOfFreedom dofw1 = eel.infoRight->getOppVertexDof(eel.dofEdge);

  WorldVector<double> evec = eel.infoLeft->getEdge(eel.dofEdge);
  WorldVector<double> hvec = eel.infoLeft->getCoordFromGlobalIndex(dofw2)
                            -eel.infoRight->getCoordFromGlobalIndex(dofw1);

  double ee = evec * evec;
  double hh = hvec * hvec;
  double eh = evec * hvec;

  double sqg= std::sqrt(ee*hh - eh*eh); // sqrt(|g|)
  // scale
  double ce = f * eh / sqg;
  double ch = -f * ee / sqg;

  rowMapper[dofv1] = -ce;
  rowMapper[dofv2] =  ce;
  rowMapper[dofw1] = -ch;
  rowMapper[dofw2] =  ch;

  return rowMapper;
}

// < f alpha , edge > = 0.5 <alpha, edge> ( <f,v1> + <f,v2> )
edgeRowValMapper AverageVertexAndEdgeVecAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  double f = 0.5 * (*evec)[eel] * fac * factor;

  rowMapper[eel.dofEdge.second] = f;
  rowMapper[eel.dofEdge.first] = f;

  return rowMapper;
}
