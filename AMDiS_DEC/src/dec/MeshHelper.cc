#include "MeshHelper.h"
#include "WorldVectorHelper.h"
#include "elVolumesInfo2d.h"

namespace AMDiS{

DOFVector<int> getConnections(const FiniteElemSpace *feSpace) {
  DOFVector<int> n(feSpace, "NrOfConnections");
  n = 0;

  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL); el; el = stack.traverseNext(el)) {
    for (int i = 0; i < 3; i++) {
      DegreeOfFreedom dof = el->getElement()->getDof(i,0);
      n[dof]++;
    }
  }

  return n;
}

DOFVector<double> get1RingVols(const FiniteElemSpace *feSpace) {
  DOFVector<double> V(feSpace,"1RingVolumes");
  V = 0.0;
  
  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_DET); el; el = stack.traverseNext(el)) {
    for (int i = 0; i < 3; i++) {
      DegreeOfFreedom dof = el->getElement()->getDof(i,0);
      V[dof] += 0.5 * el->getDet();
    }
  }

  return V;
}

DOFVector<double> getDualVols(const FiniteElemSpace *feSpace) {
  DOFVector<double> V(feSpace,"DualVolumes");
  V = 0.0;
  
  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET); el; el = stack.traverseNext(el)) {
    ElVolumesInfo2d vols(el);
    for (int i = 0; i < 3; i++) {
      DegreeOfFreedom dof = el->getElement()->getDof(i,0);
      V[dof] += vols.getDualVertexVol(i);
    }
  }

  return V;
}


DOFVector<double> getVoronoiRadii(const FiniteElemSpace *feSpace) {
  DOFVector<double> Radii(feSpace,"VoronoiRadii");
  Radii = 0.0;

  DOFVector<double> dualVols = getDualVols(feSpace);
  
  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET); el; el = stack.traverseNext(el)) {
    ElVolumesInfo2d vols(el);
    //set element radii
    ElementVector elR(3);
    elR = 0.0;
    for (int i = 0; i < 3; i++) {
       double edgeLen = vols.getOppEdgeLen((i+2)%3);
       elR[i] += edgeLen;
       elR[(i+1)%3] += edgeLen;
    }
    elR *= 0.25;
    
    for (int i = 0; i < 3; i++) {
      DegreeOfFreedom dof = el->getElement()->getDof(i,0);
      Radii[dof] += vols.getDualVertexVol(i) * elR[i] / dualVols[dof];
    }
  }

  return Radii;
}

}
