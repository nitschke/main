#include "DofVertexVector.h"
#include "DofEdgeVector.h"

using namespace AMDiS;
using namespace dec;

DofEdgeVector* DofVertexVector::exteriorDerivative() const {
  DofEdgeVector *exd = new DofEdgeVector(emesh, "d(" + getName() + ")");
  exd->set(0.0);
  
  vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
  for (; edgeIter != emesh->getEdges()->end(); ++edgeIter) {
    (*exd)[*edgeIter] = (*this)[edgeIter->dofEdge.second] - (*this)[edgeIter->dofEdge.first];
  }

  return exd;
}

DofVertexVector* DofVertexVector::laplace() const {
  DofVertexVector *lb = new DofVertexVector(emesh, "Laplace_Of_" + getName());
  lb->set(0.0);
  
  // voronoi areas on vertices
  DOFVector<double> vvol(emesh->getFeSpace(), "vvol");
  vvol.set(0.0);

  vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
  for (; edgeIter !=  emesh->getEdges()->end(); ++edgeIter) {
    double lenP = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    double lenD = edgeIter->infoLeft->getDualEdgeLen(edgeIter->dofEdge)
                + edgeIter->infoRight->getDualEdgeLen(edgeIter->dofEdge);
    DegreeOfFreedom dof1 = edgeIter->dofEdge.first;
    DegreeOfFreedom dof2 = edgeIter->dofEdge.second;
    double edgeval = (lenD / lenP) * ( (*this)[dof2] - (*this)[dof1] );
    (*lb)[dof1] += edgeval;
    (*lb)[dof2] -= edgeval;

    // this ensure that we not doubled the local voronoi areas
    vvol[dof1] += edgeIter->infoLeft->getDualVertexVol_global(dof1);
    vvol[dof2] += edgeIter->infoRight->getDualVertexVol_global(dof2);
  }

  // scaling with precalculated voronoi areas
  DOFVector<double>::Iterator lbIter(const_cast<DofVertexVector*>(lb), USED_DOFS);
  DOFVector<double>::Iterator vvolIter(const_cast<DOFVector<double>*>(&vvol), USED_DOFS);
  for (lbIter.reset(), vvolIter.reset(); !lbIter.end(); ++lbIter, ++vvolIter) {
    (*lbIter) /= (*vvolIter);
  }

  return lb;
}

DofEdgeVector* DofVertexVector::rotOnEdges_evalOnOppositeVertices() const {
  DofEdgeVector *rot = new DofEdgeVector(emesh, "rot");
  rot->set(0.0);

  vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
  for (; edgeIter != emesh->getEdges()->end(); ++edgeIter) {
    double volL = edgeIter->infoLeft->getVol();
    double volR = edgeIter->infoRight->getVol();
    double len = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    DegreeOfFreedom dofL = edgeIter->infoLeft->getOppVertexDof(edgeIter->dofEdge);
    DegreeOfFreedom dofR = edgeIter->infoRight->getOppVertexDof(edgeIter->dofEdge);
    (*rot)[*edgeIter] = - 0.5 * len * len * ((*this)[dofL] - (*this)[dofR]) / (volL + volR);
  }

  return rot;
}



DofEdgeVector* DofVertexVector::rotOnEdges_evalOnAllVertices() const {
  DofEdgeVector *rot = new DofEdgeVector(emesh, "rot");
  rot->set(0.0);

  vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
  for (; edgeIter != emesh->getEdges()->end(); ++edgeIter) {
    //double volL = edgeIter->infoLeft->getVol();
    //double volR = edgeIter->infoRight->getVol();
    //double len = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    
    // edge e = [v1,v2]
    WorldVector<double> v1v2 = edgeIter->infoLeft->getEdge(edgeIter->dofEdge);
    DegreeOfFreedom dofv1 = edgeIter->dofEdge.first;
    DegreeOfFreedom dofv2 = edgeIter->dofEdge.second;
    double dfe = (*this)[dofv2] - (*this)[dofv1];

    // "dual" edge h = [w1,w2]
    DegreeOfFreedom dofw2 = edgeIter->infoLeft->getOppVertexDof(edgeIter->dofEdge);
    DegreeOfFreedom dofw1 = edgeIter->infoRight->getOppVertexDof(edgeIter->dofEdge);
    WorldVector<double> w2 = edgeIter->infoLeft->getCoordFromGlobalIndex(dofw2);
    WorldVector<double> w1 = edgeIter->infoRight->getCoordFromGlobalIndex(dofw1);
    WorldVector<double> w1w2 = w2 - w1;
    double dfh = (*this)[dofw2] - (*this)[dofw1];

    // metric quants
    double ee = v1v2*v1v2;
    double hh = w1w2*w1w2;
    double eh = v1v2*w1w2;
    double sqg= std::sqrt(ee*hh - eh*eh); // sqrt(|g|)
    //double sqg= 2.0 * (edgeIter->infoLeft->getVol() + edgeIter->infoRight->getVol()); // ~sqrt(|g|); "=" in the flat case; slightly worse on sphere test with Rot(z);


    (*rot)[*edgeIter] = (eh*dfe - ee*dfh) / sqg;
  }

  return rot;
}

DofEdgeVector* DofVertexVector::rotOnEdges_evalOnAllVerticesAlt() const {
  DofEdgeVector *rot = new DofEdgeVector(emesh, "rot");
  rot->set(0.0);

  vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
  for (; edgeIter != emesh->getEdges()->end(); ++edgeIter) {
    //double volL = edgeIter->infoLeft->getVol();
    //double volR = edgeIter->infoRight->getVol();
    //double len = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    
    // edge e = [v1,v2]
    WorldVector<double> v1v2 = edgeIter->infoLeft->getEdge(edgeIter->dofEdge);
    DegreeOfFreedom dofv1 = edgeIter->dofEdge.first;
    DegreeOfFreedom dofv2 = edgeIter->dofEdge.second;
    double fv1 = (*this)[dofv1];
    double fv2 = (*this)[dofv2];
    double dfe = fv2 - fv1;
    double fc = 0.5 * (fv1 + fv2);

    // "dual" edges hi = s[wi,c]
    DegreeOfFreedom dofw2 = edgeIter->infoLeft->getOppVertexDof(edgeIter->dofEdge);
    DegreeOfFreedom dofw1 = edgeIter->infoRight->getOppVertexDof(edgeIter->dofEdge);
    double fw1 = (*this)[dofw1];
    double fw2 = (*this)[dofw2];
    WorldVector<double> w2 = edgeIter->infoLeft->getCoordFromGlobalIndex(dofw2);
    WorldVector<double> w1 = edgeIter->infoRight->getCoordFromGlobalIndex(dofw1);
    WorldVector<double> c = edgeIter->infoLeft->getEdgeCenter(edgeIter->dofEdge);
    WorldVector<double> cw2 = w2 - c;
    WorldVector<double> w1c = c - w1;
    double dfh1 = fc - fw1;
    double dfh2 = fw2 - fc;

    // metric quants
    double ee = v1v2*v1v2;
    double hh1 = w1c*w1c;
    double hh2 = cw2*cw2;
    double eh1 = v1v2*w1c;
    double eh2 = v1v2*cw2;
    double sqg1= std::sqrt(ee*hh1 - eh1*eh1); // sqrt(|g|)
    double sqg2= std::sqrt(ee*hh2 - eh2*eh2); // sqrt(|g|)


    (*rot)[*edgeIter] = 0.5*((eh1*dfe - ee*dfh1) / sqg1 +  (eh2*dfe - ee*dfh2) / sqg2);
  }

  return rot;
}
