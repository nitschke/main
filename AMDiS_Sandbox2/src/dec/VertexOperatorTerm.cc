#include "VertexOperatorTerm.h"
#include "ElVolumesInfo2d.h"

vertexRowValMapper LaplaceBeltramiAtVertices::evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {
  vertexRowValMapper rowMapper;

  DegreeOfFreedom vdof = (pos == FIRSTVERTEX) ? eel.dofEdge.first : eel.dofEdge.second;

  // calc voronoi vol
  double vvol = 0.0;
  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    vvol += eIter.getFace()->getDualVertexVol_global(vdof);
  }
  // factor scale
  vvol *= fac * factor;

  //build mapper
  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    //  |*e| / |e|
    double pdratio =  (  eIter->infoLeft()->getDualEdgeLen(eIter->dofEdge) 
                      + eIter->infoRight()->getDualEdgeLen(eIter->dofEdge) )
                    / eIter->infoLeft()->getEdgeLen(eIter->dofEdge);
    double c = vvol * pdratio;
    DegreeOfFreedom wdof = (eIter.pointOutward()) ? eIter->dofEdge.second : eIter->dofEdge.first;
    rowMapper[wdof] = c;
    rowMapper[vdof] = -c;
  }

  return rowMapper;
}
