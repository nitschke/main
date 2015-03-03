#include "DofEdgeVector.h"

// t in [-1,1]
inline WorldVector<double> linTrans(double t, WorldVector<double> p,  WorldVector<double> q) {
  return 0.5 * ((1.0 - t) * p + (1.0 + t) * q);
}

DofEdgeVector::interpolGL4(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
            AbstractFunction<WorldVector<double>, WorldVector<double> > *proj,
            AbstractFunction<WorldMatrix<double>, WorldVector<double> > *jproj) {

  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(&coords);

  double t1 = -sqrt(3./7. + 2./7. * sqrt(6./5.));
  double t2 = -sqrt(3./7. - 2./7. * sqrt(6./5.));
  double t3 = -t2;
  double t4 = -t1;

  double a1 = (18. - sqrt(30.)) / 36.;
  double a2 = (18. + sqrt(30.)) / 36.;
  double a3 = a2;
  double a4 = a1;

  vector<DofEdge>::iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->first];
    WorldVector<double> q = coords[edgeIter->second];
  }

}
