#ifndef DOFEDGEVECTOR_H
#define DOFEDGEVECTOR_H

#include "Dec_fwd.h"
#include "EdgeMesh.h"
#include "io/VtkVectorWriter.h"


using namespace std;

namespace AMDiS { namespace dec {


class DofEdgeVector {

public:

  DofEdgeVector() {}

  DofEdgeVector(const EdgeMesh *edgeMesh_, std::string name_) : 
      edgeMesh(edgeMesh_), name(name_), edgeVals(edgeMesh->getNumberOfEdges()) {}

  //copy-constructor (shallow mesh copy)
  DofEdgeVector(const DofEdgeVector &dev) : edgeMesh(dev.getEdgeMesh()),
                                            name(dev.getName()),
                                            edgeVals(*dev.getEdgeVector()) {}

  const vector<double>* getEdgeVector() const {return &edgeVals;}

  const EdgeMesh* getEdgeMesh() const {return edgeMesh;}

  const string getName() const {return name;}

  void setName(std::string name_) {name = name_;}

  // set alpha_d([p,q]) = <alpha, [p,q]>
  void set(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d);

  // set alpha_d(e) = <alpha, e>
  void set(AbstractFunction<double, EdgeElement > *alpha_d);

  // set -|e|/|*e| <alpha, *[p,q]> = <*alpha, [p,q]>
  //TODO: rename is more interpolating (with approx |pi(e)| and |pi(*e)|)
  void setDual(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d);

  // set  <alpha, sigma^1> = val (for all sigma^1)
  void set(double val);


  // Gauss-Legendre-Intigration (n=4) of I_edge[alpha]
  void interpolGL4(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
            AbstractFunction<WorldVector<double>, WorldVector<double> > *proj,
            AbstractFunction<WorldMatrix<double>, WorldVector<double> > *jproj);

  // lin approx of the abstract edge (= edge itself) and trapz-rule
  void interpolLinTrapz(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha);
  
  // use the midpoint rule and linear edge and assume that is a eligible interpoint, where the edge lie in the tangential space (see mean value theorem)
  void interpolLinMidpoint(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha);

  // use the midpoint rule and assume that is a eligible interpoint, where the edge lie in the tangential space (see mean value theorem)
  void interpolMidpoint(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
                        AbstractFunction<WorldVector<double>, WorldVector<double> > *proj);

  // Newton-Cotes 
  void interpolNC(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
                  int n,
                        AbstractFunction<WorldVector<double>, WorldVector<double> > *proj=NULL);

  void interpol(AbstractFunction<WorldVector<double>, WorldVector<double> > *vec);

  DOFVector< WorldVector<double> > getSharpEdgeRingLinMod();

  DOFVector< WorldVector<double> > getSharpHirani();

  DOFVector< WorldVector<double> > getSharpFaceAverage();

  // sharp procedure is set in the init-file or FaceAveraging 
  DOFVector< WorldVector<double> > getSharp(ProblemStat *ps);

  // not stable
  DOFVector< WorldVector<double> > getSharpEdgeRing();

  map<int, std::vector<double> > getSharpOnFaces();

  // - delta * d = rot*rot
  DofEdgeVector laplaceBeltrami();

  // - d * delta = grad * div
  DofEdgeVector laplaceCoBeltrami();

  DOFVector<double> divergence() const;

  // l2-Norm -> L2-Norm on K^(1)
  double l2Norm() {
    vector<double>::iterator valIter = edgeVals.begin();
    double norm2 = 0.0;
    for (; valIter != edgeVals.end(); ++valIter){
      norm2 += (*valIter) * (*valIter);
    }
    return sqrt(norm2);
  }

  // Error: l2Norm(alpha - solution) / Vol(K^(1))
  double error(const DofEdgeVector &sol) {
    DofEdgeVector errVec(*this);
    errVec -= sol;
    return errVec.l2Norm() / edgeMesh->getVol();
  }

  double absMax() {
    double max = -1.0;
    vector<double>::iterator valIter = edgeVals.begin();
    for (; valIter != edgeVals.end(); ++valIter) {
      double absVal = abs(*valIter);
      if (absVal > max) max = absVal;
    }
    return max;
  }

  double errorMax(const DofEdgeVector &sol) {
    DofEdgeVector errVec(*this);
    errVec -= sol;
    return errVec.absMax();
  }

  // rise indices resp. to local metric gPD = |e|^2 * delta_ij, 
  // i.e. <alpha,e> -> < alpha/|e|^2 , e > = <alpha^#PD,e>
  DofEdgeVector getLocalPDSharp() {
    DofEdgeVector sharpVec(edgeMesh, name + "PDSharp");
    vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
    for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter){
      WorldVector<double> eP = edgeIter->infoLeft->getEdge(edgeIter->dofEdge);
      double lenP2 = eP*eP;
      sharpVec[edgeIter->edgeDof] = edgeVals[edgeIter->edgeDof] / lenP2;
    }
    return sharpVec;
  }
  
  // be careful with the meaning of the result.
  // we approx. the integral of the vals with a e *e decomposition of the surface
  // this makes sence for e.q. scalar values, like norms, on the Edges
  // i.e. Int_M(f^2) approx Sum_edges( 0.5 * |e| * |*e| * (f_e)^2 )
  double L2NormSquared() {
    double norm2 = 0.0;
    vector<double>::iterator valIter = edgeVals.begin();
    vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
    for (; valIter != edgeVals.end(); ++valIter, ++edgeIter){
      DofEdge dofe = edgeIter->dofEdge;
      // 
      norm2 += 0.5 * edgeIter->infoLeft->getEdgeLen(dofe)  
                   * (edgeIter->infoLeft->getDualEdgeLen(dofe) + edgeIter->infoRight->getDualEdgeLen(dofe))
                   * (*valIter) * (*valIter);
    }
    return norm2;
  }

  // be careful with the meaning of the result.
  // we approx. the integral of the vals with a e *e decomposition of the surface
  // this makes sence for e.q. scalar values, like norms, on the Edges
  // i.e. Int_M(f) approx Sum_edges( 0.5 * |e| * |*e| * (f_e)^2 )
  double surfaceIntegration() {
    double norm2 = 0.0;
    vector<double>::iterator valIter = edgeVals.begin();
    vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
    for (; valIter != edgeVals.end(); ++valIter, ++edgeIter){
      DofEdge dofe = edgeIter->dofEdge;
      // 
      norm2 += 0.5 * edgeIter->infoLeft->getEdgeLen(dofe)  
                   * (edgeIter->infoLeft->getDualEdgeLen(dofe) + edgeIter->infoRight->getDualEdgeLen(dofe))
                   * (*valIter);
    }
    return norm2;
  }

  void evalFunction(AbstractFunction<double,double> *fun) {
    vector<double>::iterator valIter = edgeVals.begin();
    for (; valIter != edgeVals.end(); ++valIter) {
      *valIter = (*fun)(*valIter);
    }
  }

  DofEdgeVector& operator+=(const DofEdgeVector& a) {
    //TODO: test exits
    vector<double>::const_iterator aIter = a.getEdgeVector()->begin();
    vector<double>::iterator rvalIter = edgeVals.begin();
    for (; aIter != a.getEdgeVector()->end() || rvalIter != edgeVals.end() ; ++aIter, ++rvalIter) {
      (*rvalIter) += (*aIter);
    }
    return *this;
  }

  DofEdgeVector& operator-=(const DofEdgeVector& a) {
    //TODO: test exits
    vector<double>::const_iterator aIter = a.getEdgeVector()->begin();
    vector<double>::iterator rvalIter = edgeVals.begin();
    for (; aIter != a.getEdgeVector()->end() || rvalIter != edgeVals.end() ; ++aIter, ++rvalIter) {
      (*rvalIter) -= (*aIter);
    }
    return *this;
  }

  DofEdgeVector& operator*=(const DofEdgeVector& a) {
    //TODO: test exits
    vector<double>::const_iterator aIter = a.getEdgeVector()->begin();
    vector<double>::iterator rvalIter = edgeVals.begin();
    for (; aIter != a.getEdgeVector()->end() || rvalIter != edgeVals.end() ; ++aIter, ++rvalIter) {
      (*rvalIter) *= (*aIter);
    }
    return *this;
  }

  DofEdgeVector& operator+=(const double& c) {
    //TODO: test exits
    vector<double>::iterator rvalIter = edgeVals.begin();
    for (;rvalIter != edgeVals.end() ; ++rvalIter) {
      (*rvalIter) += c;
    }
    return *this;
  }

  DofEdgeVector& operator*=(const double& c) {
    //TODO: test exits
    vector<double>::iterator rvalIter = edgeVals.begin();
    for (;rvalIter != edgeVals.end() ; ++rvalIter) {
      (*rvalIter) *= c;
    }
    return *this;
  }

  inline double& operator[](DegreeOfFreedom dof) 
  {
    //TODO: test exits
    return edgeVals[dof];
  }

  inline double& operator[](const EdgeElement &eel) 
  {
    //TODO: test exits
    return edgeVals[eel.edgeDof];
  }

  inline double const& operator[](const EdgeElement &eel) const
  {
    //TODO: test exits
    return edgeVals[eel.edgeDof];
  }

  inline double& operator[](EdgeElement *eel) 
  {
    //TODO: test exits
    return edgeVals[eel->edgeDof];
  }


  void writeFile(string name) const;

  void writeSharpFile(string name, ProblemStat *ps) {
    DOFVector< WorldVector<double> > sharp = getSharp(ps);
    io::VtkVectorWriter::writeFile(sharp, name);
  }

  ~DofEdgeVector() {};

  
  
 
protected:

  void set(const DenseVector &vec) {
    using namespace mtl;
    int n = size(vec);
    TEST_EXIT(n == edgeVals.size())("uncool!\n");
    for (int i = 0; i < n; ++i) {
      edgeVals[i] = vec[i];
    }
  }


  const EdgeMesh *edgeMesh;

  //TODO: change to mtl-vector -> no copy operation (see private set)
  vector<double> edgeVals;

  std::string name;

  friend class DecProblemStat;
};

// holds the hodge dual *alpha in addition on the primal edge mesh; i.e. (P)rimal-(D)ual-Vals are saved
class DofEdgeVectorPD : public DofEdgeVector {
public:
  
  DofEdgeVectorPD(const EdgeMesh *edgeMesh_, std::string name_) : 
      DofEdgeVector(edgeMesh_, name_), edgeDualVals(edgeMesh->getNumberOfEdges()), isP(true) {}

  DofEdgeVectorPD(const DofEdgeVector &primal, const DofEdgeVector &dual) 
      : DofEdgeVector(primal), edgeDualVals(*dual.getEdgeVector()), isP(true) {}

  bool isPrimal() {
    return isP;
  }

  // swap(Primal Vals, Dual Vals); use std:swap
  void swapPD() {
    swap(edgeVals, edgeDualVals);
    isP = !isP;
  }

  // then edgeVals holds the primal vals
  void makePrimal() {
    if (!isP) swapPD();
  }

  // then edgeVals holds the dual vals
  void makeDual() {
    if (isP) swapPD();
  }

  //Dual is now primal
  void bakeDual() {
    swap(edgeVals, edgeDualVals);
  }

  // TODO: andere datenstruktur um kopien zu vermeiden
  DofEdgeVector getDual() {
    makeDual();
    DofEdgeVector duals(*this);
    makePrimal();
    return duals;
  }

  void set(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d,
           BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alphaDual_d) {
    makePrimal();
    DofEdgeVector::set(alpha_d);
    makeDual();
    DofEdgeVector::set(alphaDual_d);
    makePrimal();
  }

  // eval on edge on dual edge [<alpha, e> , -|e|/|*e| <alpha, *e>]
  void set(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d) {
    makePrimal();
    DofEdgeVector::set(alpha_d);
    makeDual();
    DofEdgeVector::setDual(alpha_d);
    makePrimal();
  }

  void interpol(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha);
  void interpol(AbstractFunction<WorldVector<double>, WorldVector<double> > *vec);

  void normalize(double eps = 0.0);

  DofEdgeVector getNormOnEdges() const;

  //TODO: implement DofEdgeVector with template argument 
  //      and use DofEdgeVector<WorldVector<double> > for rised indizes on Edges and file writing 
  void writeSharpOnEdgesFile(string name);

  double getDirichletEnergy(double divFac = 1.0, double rotFac = 1.0) {
    makePrimal();
    DOFVector<double> div = DofEdgeVector::divergence();
    makeDual();
    // in fact "-rot", but this is not neccesary here, because we use the squares; 
    // follows from the identity : <div *alpha, v> = -<rot alpha, v>
    DOFVector<double> rot = DofEdgeVector::divergence();
    makePrimal();

    return divFac * div.L2NormSquare() + rotFac * rot.L2NormSquare();
  }

  

private:

  vector<WorldVector<double> > getSharpVecOnEdges();

  //TODO: change to mtl-vector -> no copy operation (see private set)
  vector<double> edgeDualVals; // *alpha

  bool isP; //true if edgeVals holds the original alpha, false else
  
};


inline DofEdgeVector operator+(DofEdgeVector a, const DofEdgeVector& b) {
  return a += b;
}

inline DofEdgeVector operator-(DofEdgeVector a, const DofEdgeVector& b) {
  return a -= b;
}

inline DofEdgeVector operator*(DofEdgeVector a, const DofEdgeVector& b) {
  return a *= b;
}

inline DofEdgeVector operator*(const double& c, DofEdgeVector a) {
  return a *= c;
}

}}
#endif
