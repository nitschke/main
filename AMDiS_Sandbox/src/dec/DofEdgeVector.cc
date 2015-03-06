#include "DofEdgeVector.h"

// t in [-1,1]
inline WorldVector<double> linTrans(double t, WorldVector<double> p,  WorldVector<double> q) {
  return 0.5 * ((1.0 - t) * p + (1.0 + t) * q);
}

void DofEdgeVector::interpolGL4(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
            AbstractFunction<WorldVector<double>, WorldVector<double> > *proj,
            AbstractFunction<WorldMatrix<double>, WorldVector<double> > *jproj) {

  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<double> t(4);
  t[0] = -sqrt(3./7. + 2./7. * sqrt(6./5.));
  t[1] = -sqrt(3./7. - 2./7. * sqrt(6./5.));
  t[2] = -t[1];
  t[3] = -t[0];

  vector<double> a(4);
  a[0] = (18. - sqrt(30.)) / 36.;
  a[1] = (18. + sqrt(30.)) / 36.;
  a[2] = a[1];
  a[3] = a[0];

  vector<DofEdge>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    (*valIter) = 0.0;
    WorldVector<double> p = coords[edgeIter->first];
    WorldVector<double> q = coords[edgeIter->second];
    for (int i= 0; i < 4; i++) {
        WorldVector<double> linS = linTrans(t[i], p, q);
        // alpha_(s(t))(dt(s(t))) = alpha_(pi(s*))(J_pi(s*(t)).dt(s*(t)))
       (*valIter) += a[i] * (*alpha)((*proj)(linS), 0.5 * (*jproj)(linS) * (q - p)); 
    }
    //cout << p << " , " << q << " : " << (*valIter) << endl;
  }
}

void DofEdgeVector::interpolLinTrapz(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha) {
  
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<DofEdge>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->first];
    WorldVector<double> q = coords[edgeIter->second];
    WorldVector<double> sigma = q - p;
    (*valIter) = 0.5 * ((*alpha)(p, sigma) + (*alpha)(q, sigma));
  }
}

void DofEdgeVector::interpolLinMidpoint(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<DofEdge>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->first];
    WorldVector<double> q = coords[edgeIter->second];
    (*valIter) = (*alpha)(0.5 * (p + q), q - p);
  }
}

void DofEdgeVector::interpolMidpoint(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
                                     AbstractFunction<WorldVector<double>, WorldVector<double> > *proj) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<DofEdge>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->first];
    WorldVector<double> q = coords[edgeIter->second];
    (*valIter) = (*alpha)((*proj)(0.5 * (p + q)), q - p);
  }
}

void DofEdgeVector::interpolNC(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha, int n,
                                     AbstractFunction<WorldVector<double>, WorldVector<double> > *proj) {
  vector<double> t(n);
  vector<double> a(n);
  double h;
  
  switch(n) {
    case 3: //Simpson 
        t[0] = -1.;
        t[1] = 0.;
        t[2] = 1.;
        a[0] = 1./6.;
        a[1] = 4./6.;
        a[2] = 1./6.;
        h = 1.0;
        break;
    case 7: //Weddle
        t[0] = -1.0;
        t[1] = -2.0 / 3.0;
        t[2] = -1.0 / 3.0;
        t[3] = 0.0;
        t[4] = -t[2];
        t[5] = -t[1];
        t[6] = -t[0];
        a[0] = 41. / 840.;
        a[1] = 216. / 840.;
        a[2] = 27. / 840.;
        a[3] = 272. / 840.;
        a[4] = a[2];
        a[5] = a[1];
        a[6] = a[0];
        h = 1. /3.;
        break;
    default:
        MSG("DofEdgeVector::interpolNC not implemented for this parameter");
        return;
  }

  
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<DofEdge>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->first];
    WorldVector<double> q = coords[edgeIter->second];
    (*valIter) = 0.0;
    if (!proj) {
      for (int i = 0; i < n; i++) {
        (*valIter) += a[i] * (*alpha)(linTrans(t[i], p, q), q - p);
      }
    } else {
      vector<WorldVector<double> > points(n);
      points[0] = p;
      points[n-1] = q;
      for (int i = 1; i < n-1; i++) {
        points[i] = (*proj)(linTrans(t[i], p, q));
      }
      vector<WorldVector<double> > derivations(n);
      derivations[0] = (points[1] - p) / h; //forward
      derivations[n-1] = (q - points[n-2]) / h; //backward
      for (int i = 1; i < n-1; i++) {
        derivations[i] = 0.5 * (points[i+1] - points[i-1]) / h; //center
      }
      for (int i = 0; i < n; i++) {
        (*valIter) += 2.0 * a[i] * (*alpha)(points[i], derivations[i]);
      }
    }  
  }
}


void DofEdgeVector::set(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<DofEdge>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter) {
    (*valIter) = (*alpha_d)(coords[edgeIter->first], coords[edgeIter->second]);
  }
}

