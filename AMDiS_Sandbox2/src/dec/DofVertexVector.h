#ifndef DOFVERTEXVECTOR_H
#define DOFVERTEXVECTOR_H


#include "Dec_fwd.h"
#include "EdgeMesh.h"
#include "io/VtkVectorWriter.h"

using namespace std;

namespace AMDiS { namespace dec {

//Nur Behelfsloesung
class DofVertexVector : public DOFVector<double> {
public:

  //DofVertexVector() {}

  DofVertexVector(const EdgeMesh *edgeMesh_, std::string name_) :
      DOFVector<double>(edgeMesh_->getFeSpace(), name_), emesh(edgeMesh_) {}

  //copy-constructor (shallow mesh copy)
  DofVertexVector(const DofVertexVector &dev) : DOFVector<double>(dev), emesh(dev.getEdgeMesh()) {}

  const EdgeMesh* getEdgeMesh() const {return emesh;}

  DofEdgeVector* rotOnEdges() const {
    //return rotOnEdges_evalOnOppositeVertices();
    //return rotOnEdges_evalOnAllVertices();
    return rotOnEdges_evalOnAllVertices();
  }

  // nicht konsistent (test: sphere rot(z))
  //<Rot(f),e> = -(|e|^2 / 2(|T1|+|T2|)) (f(w2) - f(w1))
  //w1,w2...opposite vertices resp. e
  DofEdgeVector* rotOnEdges_evalOnOppositeVertices() const;


  DofEdgeVector* rotOnEdges_evalOnAllVertices() const;

  // slightly better, but more expensive
  DofEdgeVector* rotOnEdges_evalOnAllVerticesAlt() const;

  double errorMax(const DofVertexVector &sol) {
    DofVertexVector errVec(*this);
    errVec -= sol;
    return errVec.absMax();
  }

  double errorMaxRel(const DofVertexVector &sol) {
    DofVertexVector errVec(*this);
    errVec -= sol;
    return errVec.absMax() / sol.absMax();
  }

  double errorL2(const DofVertexVector &sol) {
    DofVertexVector errVec(*this);
    errVec -= sol;
    return errVec.L2NormSquare();
  }

  double errorL2Rel(const DofVertexVector &sol) {
    DofVertexVector errVec(*this);
    errVec -= sol;
    return errVec.L2NormSquare() / sol.L2NormSquare();
  }


  DofVertexVector& operator-=(const DofVertexVector& a) {
    //DOFVector<double>::Iterator xIterator(dynamic_cast<DOFVector<double>*>(this), USED_DOFS);
    //DOFVector<double>::Iterator yIterator(const_cast<DOFVector<double>*>(&a), USED_DOFS);
    DOFVector<double>::Iterator xIterator(this, USED_DOFS);
    DOFVector<double>::Iterator yIterator(const_cast<DofVertexVector*>(&a), USED_DOFS);
    
    for (xIterator.reset(), yIterator.reset(); !xIterator.end();
	              ++xIterator, ++yIterator) {
      (*xIterator) -= (*yIterator); 
    }

    return *this;
  }

private:

  const EdgeMesh *emesh;

};

inline DofVertexVector operator-(DofVertexVector a, const DofVertexVector& b) {
  return a -= b;
}

}}

#endif
