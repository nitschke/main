#include "VertexOperatorTerm.h"
#include "ElVolumesInfo2d.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"

using namespace AMDiS;
using namespace dec;


vertexRowValMapper LaplaceBeltramiAtVertices::evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {
  vertexRowValMapper rowMapper;

  DegreeOfFreedom vdof = (pos == FIRSTVERTEX) ? eel.dofEdge.first : eel.dofEdge.second;
  //TODO: default constructor for double -> 0.0 ???
  rowMapper[vdof] = 0.0;

  // calc voronoi vol
  double vvol = 0.0;
  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    vvol += eIter.getFace()->getDualVertexVol_global(vdof);
  }

  //build mapper
  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    //  |*e| / |e|
    double pdratio =  (  eIter->infoLeft->getDualEdgeLen(eIter->dofEdge) 
                      + eIter->infoRight->getDualEdgeLen(eIter->dofEdge) )
                    / eIter->infoLeft->getEdgeLen(eIter->dofEdge);
    double c = fac * factor * pdratio / vvol;
    DegreeOfFreedom wdof = (eIter.pointOutward()) ? eIter->dofEdge.second : eIter->dofEdge.first;
    rowMapper[wdof] = c;
    rowMapper[vdof] -= c;
  }

  return rowMapper;
}

vertexRowValMapper VertexVecLaplaceBeltramiAtVertices::evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {
  vertexRowValMapper rowMapper = LaplaceBeltramiAtVertices::evalRow(eel, pos, factor);
  
  DegreeOfFreedom vdof = (pos == FIRSTVERTEX) ? eel.dofEdge.first : eel.dofEdge.second;

  for (vertexRowValMapper::iterator mIter = rowMapper.begin(); mIter != rowMapper.end(); ++mIter) {
    if (func) {
      mIter->second *= (*func)((*dv)[vdof]);
    } else {
      mIter->second *= (*dv)[vdof];
    }
  }

  return rowMapper;
}

vertexRowValMapper IdentityAtVertices::evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {
  vertexRowValMapper rowMapper;

  DegreeOfFreedom vdof = (pos == FIRSTVERTEX) ? eel.dofEdge.first : eel.dofEdge.second;

  rowMapper[vdof] = fac * factor;
  
  return rowMapper;
}

vertexRowValMapper VertexVecAtVertices::evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {
  vertexRowValMapper rowMapper;

  DegreeOfFreedom vdof = (pos == FIRSTVERTEX) ? eel.dofEdge.first : eel.dofEdge.second;

  rowMapper[vdof] = fac * factor * (*dv)[vdof];
  
  return rowMapper;
}

vertexRowValMapper DivAtVertices::evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {
  vertexRowValMapper rowMapper;

  DegreeOfFreedom vdof = (pos == FIRSTVERTEX) ? eel.dofEdge.first : eel.dofEdge.second;
  //TODO: default constructor for double -> 0.0 ???
  rowMapper[eel.edgeDof] = 0.0;

  // calc voronoi vol
  double vvol = 0.0;
  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    vvol += eIter.getFace()->getDualVertexVol_global(vdof);
  }

  //build mapper
  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    //  |*e| / |e|
    double pdratio =  (  eIter->infoLeft->getDualEdgeLen(eIter->dofEdge) 
                      + eIter->infoRight->getDualEdgeLen(eIter->dofEdge) )
                    / eIter->infoLeft->getEdgeLen(eIter->dofEdge);
    double c = fac * factor * pdratio / vvol;
    rowMapper[eIter->edgeDof] =  (eIter.pointOutward()) ? c : -c;
  }

  return rowMapper;
}

vertexRowValMapper InterProdPartAtVertices::evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {
  vertexRowValMapper rowMapper;

  DegreeOfFreedom vdof = (pos == FIRSTVERTEX) ? eel.dofEdge.first : eel.dofEdge.second;
  //TODO: default constructor for double -> 0.0 ???
  rowMapper[eel.edgeDof] = 0.0;

  // calc voronoi vol (times 4)
  double vvol4 = 0.0;
  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    vvol4 += eIter.getFace()->getDualVertexVol_global(vdof);
  }
  vvol4 *= 4.0;

  for(EdgeElement::EdgeRingIterator eIter(&eel, pos); !eIter.isEnd(); ++eIter) {
    //  |*e| / |e|
    double pdratio =  (  eIter->infoLeft->getDualEdgeLen(eIter->dofEdge) 
                      + eIter->infoRight->getDualEdgeLen(eIter->dofEdge) )
                     / eIter->infoLeft->getEdgeLen(eIter->dofEdge);
    double c = fac * factor * pdratio / vvol4;
    rowMapper[eIter->edgeDof] = c * (*gamma)[*eIter];
  }

  return rowMapper;
}
