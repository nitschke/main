#ifndef DOFEDGEVECTOR_H
#define DOFEDGEVECTOR_H

#include "AMDiS.h"
#include "EdgeMesh.h"

using namespace AMDiS;
using namespace std;

class DofEdgeVector {

public:

  DofEdgeVector(const EdgeMesh *edgeMesh_, std::string name_) : 
      edgeMesh(edgeMesh_), name(name_), edgeVals(edgeMesh->getNumberOfEdges()) {}

  //copy-constructor (shallow mesh copy)
  DofEdgeVector(const DofEdgeVector &dev) : edgeMesh(dev.getEdgeMesh()),
                                            name(dev.getName()),
                                            edgeVals(*dev.getEdgeVector()) {}

  const vector<double>* getEdgeVector() const {return &edgeVals;}

  const EdgeMesh* getEdgeMesh() const {return edgeMesh;}

  const string getName() const {return name;}

  // set alpha_d([p,q]) = <alpha, [p,q]>
  void set(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d);

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

  DOFVector< WorldVector<double> > getSharpEdgeRingLinMod();

  DOFVector< WorldVector<double> > getSharpHirani();

  DOFVector< WorldVector<double> > getSharpFaceAverage();

  // not stable
  DOFVector< WorldVector<double> > getSharpEdgeRing();

  map<int, std::vector<double> > getSharpOnFaces();

  // - delta * d = rot*rot
  DofEdgeVector laplaceBeltrami();

  // - d * delta = grad * div
  DofEdgeVector laplaceCoBeltrami();

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

  DofEdgeVector& operator-=(const DofEdgeVector& a) {
    //TODO: test exits
    vector<double>::const_iterator aIter = a.getEdgeVector()->begin();
    vector<double>::iterator rvalIter = edgeVals.begin();
    for (; aIter != a.getEdgeVector()->end() || rvalIter != edgeVals.end() ; ++aIter, ++rvalIter) {
      (*rvalIter) -= (*aIter);
    }
    return *this;
  }

  inline double& operator[](DegreeOfFreedom dof) 
  {
    //TODO: test exits
    return edgeVals[dof];
  }


  void writeFile(string name) const;
  
  
 
private:
  const EdgeMesh *edgeMesh;

  vector<double> edgeVals;

  std::string name;

};


DofEdgeVector operator-(DofEdgeVector a, const DofEdgeVector& b) {
  return a -= b;
}

#endif

