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

  const vector<double>* getEdgeVector() const {return &edgeVals;}

  const EdgeMesh* getEdgeMesh() const {return edgeMesh;}

  // set alpha_d([p,q]) = <alpha, [p,q]>
  void set(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d);

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



  // L2-Norm on K^(1)
  double L2Norm() {
    vector<double>::iterator valIter = edgeVals.begin();
    double norm2 = 0.0;
    for (; valIter != edgeVals.end(); ++valIter){
      norm2 += (*valIter) * (*valIter);
    }
    return sqrt(norm2);
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

  
 
private:
  const EdgeMesh *edgeMesh;

  vector<double> edgeVals;

  std::string name;

};


DofEdgeVector operator-(DofEdgeVector a, const DofEdgeVector& b) {
  return a -= b;
}

#endif
