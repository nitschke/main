#include "DofEdgeVector.h"

using namespace AMDiS;
using namespace dec;


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
  t[0] = -std::sqrt(3./7. + 2./7. * std::sqrt(6./5.));
  t[1] = -std::sqrt(3./7. - 2./7. * std::sqrt(6./5.));
  t[2] = -t[1];
  t[3] = -t[0];

  vector<double> a(4);
  a[0] = (18. - std::sqrt(30.)) / 36.;
  a[1] = (18. + std::sqrt(30.)) / 36.;
  a[2] = a[1];
  a[3] = a[0];
  
  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    (*valIter) = 0.0;
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
    for (int i= 0; i < 4; i++) {
        WorldVector<double> linS = linTrans(t[i], p, q);
        // alpha_(s(t))(dt(s(t))) = alpha_(pi(s*))(J_pi(s*(t)).dt(s*(t)))
       (*valIter) += a[i] * (*alpha)((*proj)(linS), 0.5 * (*jproj)(linS) * (q - p)); 
    }
  }
}

void DofEdgeVector::interpolLinTrapz(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha) {
  
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
    WorldVector<double> sigma = q - p;
    (*valIter) = 0.5 * ((*alpha)(p, sigma) + (*alpha)(q, sigma));
  }
}

void DofEdgeVector::interpolLinMidpoint(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
    (*valIter) = (*alpha)(0.5 * (p + q), q - p);
  }
}

void DofEdgeVector::interpolMidpoint(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
                                     AbstractFunction<WorldVector<double>, WorldVector<double> > *proj) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
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

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter){
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
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


void DofEdgeVector::interpol(AbstractFunction<WorldVector<double>, WorldVector<double> > *vec) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIterP = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIterP != edgeVals.end()  ; 
                ++edgeIter, ++valIterP) {
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
    WorldVector<double> ce = 0.5 * (p + q); // circumcenter of the primal edge
    WorldVector<double> vecval = (*vec)(ce);
    WorldVector<double> e = q - p; // primal edge = q - p
    // set primal val
    (*valIterP) = vecval * e;
  }
}


void DofEdgeVector::set(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter) {
    (*valIter) = (*alpha_d)(coords[edgeIter->dofEdge.first], coords[edgeIter->dofEdge.second]);
  }
}

void DofEdgeVector::set(AbstractFunction<double, EdgeElement > *alpha_d) {
  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter) {
    (*valIter) = (*alpha_d)(*edgeIter);
  }
}

void DofEdgeVector::setDual(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha_d) {
  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIter = edgeVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIter != edgeVals.end(); ++edgeIter, ++valIter) {
    WorldVector<double> p = edgeIter->infoRight->getCircumcenter();
    WorldVector<double> q = edgeIter->infoLeft->getCircumcenter();
    double lenP = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    double lenD = edgeIter->infoLeft->getDualEdgeLen(edgeIter->dofEdge)
                + edgeIter->infoRight->getDualEdgeLen(edgeIter->dofEdge);
    (*valIter) = - lenP * (*alpha_d)(p, q) / lenD;
  }
}

void DofEdgeVector::set(double val) {
  vector<double>::iterator valIter = edgeVals.begin();
  for (; valIter != edgeVals.end(); ++valIter) {
    (*valIter) = val;
  }
}

DOFVector< WorldVector<double> > DofEdgeVector::getSharpEdgeRingLinMod() const{
  using namespace mtl;
  DOFVector< WorldVector<double> > sharp(edgeMesh->getFeSpace(), name + "Sharp");
  DOFVector< list<EdgeElement> > edgeRings = edgeMesh->getEdgeRings();
 
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  DOFVector< WorldVector<double> >::Iterator sharpIter(&sharp, USED_DOFS);
  DOFVector< list<EdgeElement> >::Iterator ringIter(&edgeRings, USED_DOFS);
  for (sharpIter.reset(), ringIter.reset(); !sharpIter.end() ; ++sharpIter, ++ringIter) {
    dense2D<double> X(ringIter->size(), 3);
    dense_vector<double> b(ringIter->size());
    list<EdgeElement>::iterator edgeIter = ringIter->begin();
    for (int i = 0; edgeIter != ringIter->end(); ++i, ++edgeIter) {
      b[i] = edgeVals[edgeIter->edgeDof];
      WorldVector<double> sigma1 = coords[edgeIter->dofEdge.second] - coords[edgeIter->dofEdge.first];
      for (int k= 0; k < 3; k++) X[i][k] = sigma1[k];
    }
    dense2D<double> A(3,3);
    dense_vector<double> y(3), sharpwv(3);
    A = trans(X) * X;
    y = trans(X) * b;
    sharpwv = lu_solve(A, y);
    for (int k= 0; k < 3; k++) (*sharpIter)[k] = sharpwv[k];
  }
  
  return sharp; 
}

inline int getLocIndex(Element *el, DegreeOfFreedom dof) {
  for (int i = 0; i < 3; i++) if ( dof == el->getDof(i,0) ) return i;
  MSG("DofEdgeVector::getLocIndex irgendwas laeuft hier schief!");
}

DOFVector< WorldVector<double> > DofEdgeVector::getSharpHirani() {
  DOFVector< WorldVector<double> > sharp(edgeMesh->getFeSpace(), name + "Sharp");
  DOFVector< list<EdgeElement> > edgeRings = edgeMesh->getEdgeRings();
  
  WorldVector<double> faceSum;
  DOFAdmin *admin = edgeMesh->getFeSpace()->getAdmin();
  for (DegreeOfFreedom i = 0; i < sharp.getSize() ; i++) {
    if (!admin->isDofFree(i)) {
      sharp[i] = 0.0;
      list<EdgeElement>::iterator edgeIter = edgeRings[i].begin();
      for (; edgeIter != edgeRings[i].end(); ++edgeIter) {
        faceSum = 0.0;
        int locILeft = getLocIndex(edgeIter->infoLeft->getElInfo()->getElement(), i);
        faceSum += (edgeIter->infoLeft->getDualVertexVol(locILeft) / edgeIter->infoLeft->getVol()) 
                          * (edgeIter->infoLeft->getElInfo()->getGrdLambda())[locILeft];
        int locIRight = getLocIndex(edgeIter->infoRight->getElInfo()->getElement(), i);
        faceSum += (edgeIter->infoRight->getDualVertexVol(locIRight) / edgeIter->infoRight->getVol()) 
                          * (edgeIter->infoRight->getElInfo()->getGrdLambda())[locIRight];
        faceSum *= edgeVals[edgeIter->edgeDof];
        if (i == edgeIter->dofEdge.first) faceSum *= -1.0;
        sharp[i] += 0.5 * faceSum;
      }
    }
  }
  return sharp;
}

DOFVector< WorldVector<double> > DofEdgeVector::getSharpFaceAverage() const{
  using namespace mtl;
  DOFVector< WorldVector<double> > sharp(edgeMesh->getFeSpace(), name + "Sharp");
  DOFVector< list<EdgeElement> > edgeRings = edgeMesh->getEdgeRings();
 
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  //cout << *edgeMesh << endl;
  //DegreeOfFreedom coutDOF = 50401;
  //DegreeOfFreedom coutDOF = 4;
  //bool verbose = false;


  DOFAdmin *admin = edgeMesh->getFeSpace()->getAdmin();
  for (DegreeOfFreedom dof = 0; dof < sharp.getSize() ; dof++) {
    if (!admin->isDofFree(dof)) {
      //verbose = (dof == coutDOF) ? true : false;
      sharp[dof] = 0.0;
      
      WorldVector<double> v = coords[dof];
      //if (verbose) cout << "*** DOF: " << dof << " (" << v << ") ***" << endl;
      
      double voronoiVol = 0.0;
      list<EdgeElement>::iterator edgeIter = edgeRings[dof].begin();
      for (; edgeIter != edgeRings[dof].end(); ++edgeIter) {
        //if (verbose) cout << "--- Edge: " << *edgeIter;
        ElVolumesInfo2d *face = (edgeIter->dofEdge.first == dof) ? edgeIter->infoLeft : edgeIter->infoRight;
        int i = face->getLocal(dof);
        int ii = (i+1)%3;
        int iii = (i+2)%3;
        //cout << "*locLeft: " << i << " " << ii << " " << iii << endl;
        
        //find next edge with the same face in the ring TODO: improve
        list<EdgeElement>::iterator edge1it = edgeRings[dof].begin();
        while ( (edge1it->infoLeft != face && edge1it->infoRight != face) || edge1it == edgeIter ) ++edge1it;
        //if (verbose) cout << "+++ Edge: " << *edge1it << endl;


        double alpha0 = edgeVals[edgeIter->edgeDof];
        if (edgeIter->dofEdge.second == dof) alpha0 *= -1.0;
        double alpha1 = edgeVals[edge1it->edgeDof];
        if (edge1it->dofEdge.second == dof) alpha1 *= -1.0;

        WorldVector<double> edge0 = face->getElInfo()->getCoord(ii) - v;
        WorldVector<double> edge1 = face->getElInfo()->getCoord(iii) - v;

        double len02 = face->getOppEdgeLen(iii) * face->getOppEdgeLen(iii);
        double len12 = face->getOppEdgeLen(ii) * face->getOppEdgeLen(ii);
        double scal = edge0 * edge1;
        double det = len02 * len12 - scal * scal;

        double a0 = (len12 * alpha0 - scal * alpha1) / det;
        double a1 = (len02 * alpha1 - scal * alpha0) / det;


        double locVoronoiVol = face->getDualVertexVol(i); //VoronoiAverage
        //double locVoronoiVol = 1.0; //unweighted average
        //double locVoronoiVol = face->getSin(i) / (face->getOppEdgeLen(ii) * face->getOppEdgeLen(iii)); //sphere average
        
        //if (verbose) cout << (a0 * edge0 + a1 * edge1) << endl;
        sharp[dof] += locVoronoiVol * (a0 * edge0 + a1 * edge1);
        voronoiVol += locVoronoiVol;
      }
      sharp[dof] *= 1.0 / voronoiVol;
    }
  }
  return sharp; 
}


DOFVector< WorldVector<double> > DofEdgeVector::getSharpEdgeRing() {
  using namespace mtl;
  DOFVector< WorldVector<double> > sharp(edgeMesh->getFeSpace(), name + "Sharp");
  DOFVector< list<EdgeElement> > edgeRings = edgeMesh->getEdgeRings();
 
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  DOFAdmin *admin = edgeMesh->getFeSpace()->getAdmin();
  for (DegreeOfFreedom dof = 0; dof < sharp.getSize() ; dof++) {
    if (!admin->isDofFree(dof)) {
      int n = edgeRings[dof].size(); // Number of edges
      dense2D<double> G(n,n); // Gramian matrix of the edge vectors
      dense_vector<double> alpha(n); // discrete 1-form on edges

      // TODO: respect the symmetry to reduce the effort
      list<EdgeElement>::iterator edgeIter1 = edgeRings[dof].begin();
      for (int i1 = 0; edgeIter1 != edgeRings[dof].end(); ++edgeIter1, ++i1) {
        alpha[i1] = edgeVals[edgeIter1->edgeDof]; // < alpha , sigma^1_i1 >

        WorldVector<double> e1 = coords[edgeIter1->dofEdge.second] - coords[edgeIter1->dofEdge.first]; //  sigma^1_i1 (vector)
        list<EdgeElement>::iterator edgeIter2 = edgeRings[dof].begin();
        for (int i2 = 0; edgeIter2 != edgeRings[dof].end(); ++edgeIter2, ++i2) {
          WorldVector<double> e2 = coords[edgeIter2->dofEdge.second] - coords[edgeIter2->dofEdge.first]; //  sigma^1_i2 (vector)
          G[i1][i2] = e1 * e2; // (sigma^1_i1).(sigma^1_i2)
        }
      }

      dense_vector<double> a(n); // coeff. vector: alpha^# = a_i * sigma^1_i
      a = lu_solve(G,alpha);
      //cout  << G << endl;
      dense_vector<double> eig(n);
      //eig = eigenvalue_symmetric(G);
      //eig = qr_algo(G,20);
      //cout << "cond: " << (eig[0] / eig[n-1]) << endl;

      sharp[dof] = 0.0;
      list<EdgeElement>::iterator edgeIter = edgeRings[dof].begin();
      for (int i = 0; edgeIter != edgeRings[dof].end(); ++edgeIter, ++i) {
        sharp[dof] += a[i] * (coords[edgeIter->dofEdge.second] - coords[edgeIter->dofEdge.first]);
      }
    }
  }

  return sharp;
}

DOFVector< WorldVector<double> > DofEdgeVector::getSharp(ProblemStat *ps) {
  FUNCNAME("DofEdgeVector::getSharp");
  unsigned int type = 1;
  Parameters::get(ps->getName() + "->output->sharpType", type);
  TEST_EXIT(type > 0 && type < 4)("invalide sharpType");
  switch (type) {
    case 1 : return getSharpFaceAverage();
    case 2 : return getSharpEdgeRingLinMod();
    case 3 : return getSharpHirani();
  }
}

map<int, std::vector<double> > DofEdgeVector::getSharpOnFaces() {
  map<int, std::vector<double> > sharp;
  map<int, pair<ElVolumesInfo2d*, vector<EdgeElement> > > faceEdges = edgeMesh->getFaceEdges();

  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  for (map<int, pair<ElVolumesInfo2d*, vector<EdgeElement> > >::const_iterator faceEdgesIter = faceEdges.begin();
        faceEdgesIter != faceEdges.end(); ++faceEdgesIter) {
    ElVolumesInfo2d* face = (faceEdgesIter->second).first;
    WorldVector<double> elSharp;
    elSharp = 0.0;
    double weightSum = 0.0;
    for (int k = 0; k < 3; ++k){
      EdgeElement e0 = ((faceEdgesIter->second).second)[k];
      EdgeElement e1 = ((faceEdgesIter->second).second)[(k+1)%3];

      double alpha0 = edgeVals[e0.edgeDof];
      double alpha1 = edgeVals[e1.edgeDof];
    
      WorldVector<double> edge0 = coords[e0.dofEdge.second] - coords[e0.dofEdge.first];
      WorldVector<double> edge1 = coords[e1.dofEdge.second] - coords[e1.dofEdge.first];


      double len02 = edge0 * edge0;
      double len12 = edge1 * edge1;
      double scal = edge0 * edge1;
      double det = len02 * len12 - scal * scal;

      double a0 = (len12 * alpha0 - scal * alpha1) / det;
      double a1 = (len02 * alpha1 - scal * alpha0) / det;

      double locWeight = 1.0;
      //if (e0.dofEdge.first == e1.dofEdge.first || e0.dofEdge.first == e1.dofEdge.second) {
      //  locWeight = face->getDualVertexVol(face->getLocal(e0.dofEdge.first));
      //} else {
      //  locWeight = face->getDualVertexVol(face->getLocal(e0.dofEdge.second));
      //}
      weightSum += locWeight;

      elSharp += locWeight * (a0 * edge0 + a1 * edge1);
    }
    elSharp *= 1.0 / weightSum;
    //cout << elSharp << endl;
    sharp[faceEdgesIter->first] = vector<double>(elSharp.begin(), elSharp.end());
    //cout << sharp[faceEdgesIter->first][0] << " " << sharp[faceEdgesIter->first][1] << " " << sharp[faceEdgesIter->first][2] << endl;
  }
  return sharp;
}


DofEdgeVector DofEdgeVector::laplaceBeltrami() {
  DofEdgeVector lB(edgeMesh, "LaplaceBeltrami");
  //lB.set(0.0);
  
  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    // Left + (I've got the PMA)
    ElVolumesInfo2d *T = edgeIter->infoLeft;
    double dualEdgeLenLeft = T->getDualOppEdgeLen(T->getOppVertexLocal(edgeIter->dofEdge));

    double sumOnT = edgeVals[edgeIter->edgeDof];
    EdgeElement *e1 = edgeIter->edgesLeft.first;
    EdgeElement *e2 = edgeIter->edgesLeft.second;
    sumOnT += (e1->infoLeft == T) ? edgeVals[e1->edgeDof] : -edgeVals[e1->edgeDof];
    sumOnT += (e2->infoLeft == T) ? edgeVals[e2->edgeDof] : -edgeVals[e2->edgeDof];
    lB[edgeIter->edgeDof] = sumOnT / T->getVol();

    // Right -
    T = edgeIter->infoRight;
    int oppVertexLocalR = T->getOppVertexLocal(edgeIter->dofEdge);
    double dualEdgeLenRight = T->getDualOppEdgeLen(oppVertexLocalR);

    sumOnT = edgeVals[edgeIter->edgeDof];
    e1 = edgeIter->edgesRight.first;
    e2 = edgeIter->edgesRight.second;
    sumOnT += (e1->infoRight == T) ? edgeVals[e1->edgeDof] : -edgeVals[e1->edgeDof];
    sumOnT += (e2->infoRight == T) ? edgeVals[e2->edgeDof] : -edgeVals[e2->edgeDof];
    lB[edgeIter->edgeDof] += sumOnT / T->getVol();

    //scale
    lB[edgeIter->edgeDof] *= - T->getOppEdgeLen(oppVertexLocalR) / (dualEdgeLenLeft + dualEdgeLenRight);
  }

  return lB;
}

//TODO: improve (everthing is do twice)
DofEdgeVector DofEdgeVector::laplaceCoBeltrami() {
  DofEdgeVector lCB(edgeMesh, "LaplaceCoBeltrami");
  //lCB.set(0.0);
  
  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    // q-part on [p,q]
    DegreeOfFreedom vdof = edgeIter->dofEdge.second;
    double sumv = 0.0;
    double voronoiVol = 0.0;
    for (EdgeElement::EdgeRingIterator ringIter(&(*edgeIter), SECONDVERTEX); !ringIter.isEnd(); ++ringIter) {
      double l = ringIter->infoLeft->getEdgeLen(ringIter->dofEdge);
      double lStar =   ringIter->infoLeft->getDualEdgeLen(ringIter->dofEdge)
                     + ringIter->infoRight->getDualEdgeLen(ringIter->dofEdge);
      double eval = (lStar / l) * edgeVals[ringIter->edgeDof];
      sumv += (ringIter->dofEdge.first == vdof) ? eval : -eval;
      
      // TODO: not double the vol
      voronoiVol +=   ringIter->infoLeft->getDualVertexVol(ringIter->infoLeft->getLocal(vdof))
                    + ringIter->infoRight->getDualVertexVol(ringIter->infoRight->getLocal(vdof));
    }
    voronoiVol *= 0.5;
    lCB[edgeIter->edgeDof] = sumv / voronoiVol;

    // p-part on [p,q]
    vdof = edgeIter->dofEdge.first;
    sumv = 0.0;
    voronoiVol = 0.0;
    for (EdgeElement::EdgeRingIterator ringIter(&(*edgeIter), FIRSTVERTEX); !ringIter.isEnd(); ++ringIter) {
      double l = ringIter->infoLeft->getEdgeLen(ringIter->dofEdge);
      double lStar =   ringIter->infoLeft->getDualEdgeLen(ringIter->dofEdge)
                     + ringIter->infoRight->getDualEdgeLen(ringIter->dofEdge);
      double eval = (lStar / l) * edgeVals[ringIter->edgeDof];
      sumv += (ringIter->dofEdge.first == vdof) ? -eval : eval;
      
      // TODO: not double the vol
      voronoiVol +=   ringIter->infoLeft->getDualVertexVol(ringIter->infoLeft->getLocal(vdof))
                    + ringIter->infoRight->getDualVertexVol(ringIter->infoRight->getLocal(vdof));
    }
    voronoiVol *= 0.5;
    lCB[edgeIter->edgeDof] += sumv / voronoiVol;
  }

  return lCB;
}

DofEdgeVector DofEdgeVector::divOnEdgeCenter_averageOnEdgeVertices() const {
  DofEdgeVector div(edgeMesh, "div");
  div.set(0.0);

  DofVertexVector divOnVs = divergence();

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    div[*edgeIter] = 0.5 * (divOnVs[edgeIter->dofEdge.first] + divOnVs[edgeIter->dofEdge.second]);
  }

  return div;
}

DofEdgeVector DofEdgeVector::exteriorDerivativOfNorm2_oppositeEdges() const {
  using namespace mtl;

  DofEdgeVector dnorm(edgeMesh, "dnorm2");
  dnorm.set(0.0);

  std::vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    EdgeElement* e11 = edgeIter->edgesRight.first;
    EdgeElement* e12 = edgeIter->edgesLeft.first;
    EdgeElement* e21 = edgeIter->edgesRight.second;
    EdgeElement* e22 = edgeIter->edgesLeft.second;

    WorldVector<double> ev11 = e11->infoLeft->getEdge(e11->dofEdge);
    WorldVector<double> ev12 = e12->infoLeft->getEdge(e12->dofEdge);
    WorldVector<double> ev21 = e21->infoLeft->getEdge(e21->dofEdge);
    WorldVector<double> ev22 = e22->infoLeft->getEdge(e22->dofEdge);

    dense2D<double> gSecond(2,2);
    dense2D<double> gFirst(2,2);
    gFirst[0][0] = ev11 * ev11;
    gFirst[0][1] = gFirst[1][0] = ev11 * ev12;
    gFirst[1][1] = ev12 * ev12;
    gSecond[0][0] = ev21 * ev21;
    gSecond[0][1] = gSecond[1][0] = ev21 * ev22;
    gSecond[1][1] = ev22 * ev22;

    dense_vector<double> vFirst(2, 0.0);
    dense_vector<double> vSecond(2, 0.0);
    vFirst[0]  = edgeVals[e11->edgeDof];
    vFirst[1]  = edgeVals[e12->edgeDof];
    vSecond[0] = edgeVals[e21->edgeDof];
    vSecond[1] = edgeVals[e22->edgeDof];

    double normFirst  = dot((dense_vector<double>)(matrix::inv(gFirst)  * vFirst) , vFirst);
    double normSecond = dot((dense_vector<double>)(matrix::inv(gSecond) * vSecond), vSecond);

    dnorm[*edgeIter] = normSecond - normFirst;
  }

  return dnorm;
}

DofEdgeVector DofEdgeVector::exteriorDerivativOfNorm2_nextEdges() const {
  using namespace mtl;

  DofEdgeVector dnorm(edgeMesh, "dnorm2");
  dnorm.set(0.0);

  std::vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    EdgeElement* e1 = edgeIter->edgesLeft.first;
    EdgeElement* e2 = edgeIter->edgesRight.second;

    WorldVector<double> ev1 = e1->infoLeft->getEdge(e1->dofEdge);
    WorldVector<double> ev2 = e2->infoLeft->getEdge(e2->dofEdge);
    WorldVector<double> ev = edgeIter->infoLeft->getEdge(edgeIter->dofEdge);

    dense2D<double> gSecond(2,2);
    dense2D<double> gFirst(2,2);
    gFirst[0][0] = gSecond[0][0] = ev * ev;
    gFirst[0][1] = gFirst[1][0] = ev * ev1;
    gSecond[0][1] = gSecond[1][0] = ev * ev2;
    gFirst[1][1] = ev1 * ev1;
    gSecond[1][1] = ev2 * ev2;

    dense2D<double> gSecondInv = matrix::inv(gSecond);
    dense2D<double> gFirstInv = matrix::inv(gFirst);

    double det1 = gFirst[0][0] * gFirst[1][1] -  gFirst[0][1] * gFirst[1][0];
    double det2 = gSecond[0][0] * gSecond[1][1] -  gSecond[0][1] * gSecond[1][0];
    double det1Inv = 1. / det1;
    double det2Inv = 1. / det2;
    gFirstInv[0][0] = det1Inv * gFirst[1][1];
    gFirstInv[1][1] = det1Inv * gFirst[0][0];
    gFirstInv[1][0] = - det1Inv * gFirst[1][0];
    gFirstInv[0][1] = - det1Inv * gFirst[0][1];
    gSecondInv[0][0] = det2Inv * gSecond[1][1];
    gSecondInv[1][1] = det2Inv * gSecond[0][0];
    gSecondInv[1][0] = - det2Inv * gSecond[1][0];
    gSecondInv[0][1] = - det2Inv * gSecond[0][1];

    dense_vector<double> vFirst(2, 0.0);
    dense_vector<double> vSecond(2, 0.0);
    vFirst[0]  = vSecond[0] = edgeVals[edgeIter->edgeDof];
    vFirst[1]  = edgeVals[e1->edgeDof];
    vSecond[1] = edgeVals[e2->edgeDof];


    double normFirst  = dot((dense_vector<double>)(gFirstInv  * vFirst) , vFirst);
    double normSecond = dot((dense_vector<double>)(gSecondInv * vSecond), vSecond);

    //double z1 = (edgeIter->infoLeft->getCoordFromGlobalIndex(edgeIter->dofEdge.first))[2];
    //double z2 = (edgeIter->infoLeft->getCoordFromGlobalIndex(edgeIter->dofEdge.second))[2];
    //if (abs(z1) > 0.70 && abs(z1) < 0.71) {
    //cout << (1. - z1*z1) << " -> " << normFirst << endl;
    //cout << (1. - z2*z2) << " -> " << normSecond << endl << endl;
    ////WAIT;
    //}

    dnorm[*edgeIter] = normSecond - normFirst;
  }

  return dnorm;
}

DofEdgeVector DofEdgeVector::exteriorDerivativOfNorm2_FaceAverage() const {
  using namespace mtl;

  DofEdgeVector dnorm(edgeMesh, "dnorm2");
  dnorm.set(0.0);
  
  DOFVector< WorldVector<double> > vecOnVs = getSharpFaceAverage();
  //DOFVector< WorldVector<double> > vecOnVs = getSharpEdgeRingLinMod();

  std::vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    WorldVector<double> vec1 = vecOnVs[edgeIter->dofEdge.first]; 
    WorldVector<double> vec2 = vecOnVs[edgeIter->dofEdge.second]; 

    double normFirst = vec1 * vec1;
    double normSecond = vec2 * vec2;

    dnorm[*edgeIter] = normSecond - normFirst;
  }

  return dnorm;
}



DofVertexVector DofEdgeVector::divergence() const {
  DofVertexVector div(edgeMesh, name + "Div");
  div.set(0.0);
  
  // voronoi areas on vertices
  DOFVector<double> vvol(edgeMesh->getFeSpace(), "vvol");
  vvol.set(0.0);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::const_iterator valIter = edgeVals.begin();
  for (; edgeIter !=  edgeMesh->getEdges()->end(); ++edgeIter, ++valIter) {
    double lenP = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    double lenD = edgeIter->infoLeft->getDualEdgeLen(edgeIter->dofEdge)
                + edgeIter->infoRight->getDualEdgeLen(edgeIter->dofEdge);
    double divval = (lenD / lenP) * (*valIter);
    DegreeOfFreedom dof1 = edgeIter->dofEdge.first;
    DegreeOfFreedom dof2 = edgeIter->dofEdge.second;
    div[dof1] += divval;
    div[dof2] -= divval;

    // this ensure that we not doubled the local voronoi areas
    vvol[dof1] += edgeIter->infoLeft->getDualVertexVol_global(dof1);
    vvol[dof2] += edgeIter->infoRight->getDualVertexVol_global(dof2);
  }

  // scaling with precalculated voronoi areas
  DOFVector<double>::Iterator divIter(const_cast<DofVertexVector*>(&div), USED_DOFS);
  DOFVector<double>::Iterator vvolIter(const_cast<DOFVector<double>*>(&vvol), USED_DOFS);
  for (divIter.reset(), vvolIter.reset(); !divIter.end(); ++divIter, ++vvolIter) {
    (*divIter) /= (*vvolIter);
  }

  return div;
}


DofVertexVector  DofEdgeVector::averageEdgeCentersToVertices() const {
  DofVertexVector vv(edgeMesh, name);
  vv.set(0.0);

  //          __                  
  //         \      w                          
  //         /__ e>v e
  DOFVector<double> sumws(edgeMesh->getFeSpace(), "sumws");
  sumws.set(0.0);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::const_iterator valIter = edgeVals.begin();
  for (; edgeIter !=  edgeMesh->getEdges()->end(); ++edgeIter, ++valIter) {
    double len1 = 1.0 / edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    DegreeOfFreedom dof1 = edgeIter->dofEdge.first;
    DegreeOfFreedom dof2 = edgeIter->dofEdge.second;
    
    vv[dof1] += len1 * (*valIter);
    vv[dof2] += len1 * (*valIter);

    sumws[dof1] += len1;
    sumws[dof2] += len1;
  }

  //scaling
  DOFVector<double>::Iterator vvIter(const_cast<DofVertexVector*>(&vv), USED_DOFS);
  DOFVector<double>::Iterator sumwsIter(const_cast<DOFVector<double>*>(&sumws), USED_DOFS);
  for (vvIter.reset(), sumwsIter.reset(); !vvIter.end(); ++vvIter, ++sumwsIter) {
    (*vvIter) /= (*sumwsIter);
  }


  return vv;
}

DofEdgeVector DofEdgeVector::hodgeDual_unweighted() const {
  DofEdgeVector hodge(edgeMesh, "hodge dual");
  hodge.set(0.0);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    WorldVector<double> edge = edgeIter->infoLeft->getEdge(edgeIter->dofEdge);
    double edgeLen2 = edge *edge;

    double weightSum = 0.0;

    list< EdgeElement* > relatedEdges;
    relatedEdges.push_back(edgeIter->edgesLeft.first);
    relatedEdges.push_back(edgeIter->edgesLeft.second);
    relatedEdges.push_back(edgeIter->edgesRight.first);
    relatedEdges.push_back(edgeIter->edgesRight.second);

    list< double > signs;
    signs.push_back((edgeIter->dofEdge.first  == edgeIter->edgesLeft.first->dofEdge.first) ? 1.0 : -1.0);
    signs.push_back((edgeIter->dofEdge.second == edgeIter->edgesLeft.second->dofEdge.first) ? 1.0 : -1.0);
    signs.push_back((edgeIter->dofEdge.first  == edgeIter->edgesRight.first->dofEdge.second) ? 1.0 : -1.0);
    signs.push_back((edgeIter->dofEdge.second == edgeIter->edgesRight.second->dofEdge.second) ? 1.0 : -1.0);

    list< EdgeElement* >::const_iterator relEdgeIter = relatedEdges.begin();
    list< double >::const_iterator signIter = signs.begin();
    for(; relEdgeIter != relatedEdges.end(); ++relEdgeIter, ++signIter) {
      double weight = 1.0;
      weightSum += weight;

      WorldVector<double> relEdge = (*relEdgeIter)->infoLeft->getEdge((*relEdgeIter)->dofEdge);
      double relEdgeLen2 = relEdge * relEdge;
      double prodEdgeRelEdge = edge * relEdge;

      double detG = edgeLen2 * relEdgeLen2 - prodEdgeRelEdge * prodEdgeRelEdge;

      hodge[*edgeIter] += *(signIter) * weight * (prodEdgeRelEdge * edgeVals[edgeIter->edgeDof] - edgeLen2 * edgeVals[(*relEdgeIter)->edgeDof]) / std::sqrt(detG);
    }
    hodge[*edgeIter] /= weightSum;
  }
  
  return hodge;
}

void DofEdgeVector::writeFile(string name) const {
    boost::iostreams::filtering_ostream file;
    file.push(boost::iostreams::file_descriptor_sink(name, std::ios::trunc));
  file << "<?xml version=\"1.0\"?>\n";
	file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	file << "  <UnstructuredGrid>\n";
	file << "    <Piece NumberOfPoints=\"" << edgeMesh->getFeSpace()->getMesh()->getNumberOfVertices()
	     << "\" NumberOfCells=\"" << edgeVals.size() << "\">\n";
	file << "      <Points>\n";
	file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
    edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);
    DOFAdmin *admin = edgeMesh->getFeSpace()->getAdmin();
    for (DegreeOfFreedom dof = 0; dof < coords.getSize() ; dof++) {
      if (!admin->isDofFree(dof)) {
        file << coords[dof] << endl;
      }
    }
  file << "        </DataArray>\n";
	file << "      </Points>\n";
	file << "      <Cells>\n";
	file << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";
    for (int i = 1; i <= edgeVals.size(); ++i) {
      file << (i *  2) << endl;
    }
  file << "        </DataArray>\n";
	file << "        <DataArray type=\"UInt8\" Name=\"types\">\n";
    for (int i = 1; i <= edgeVals.size(); ++i) {
      file << 3 << endl;
    }
  file << "        </DataArray>\n";
	file << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";
    vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
    for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
      file << edgeIter->dofEdge.first << " " << edgeIter->dofEdge.second << endl;
    }
  file << "        </DataArray>\n";    
	file << "      </Cells>\n";
	file << "      <CellData>\n";
  file << "        <DataArray type=\"Float32\" Name=\""
	     << this->name
	     << "\" format=\"ascii\">\n";
    vector<double>::const_iterator valIter = edgeVals.begin();
    for (; valIter != edgeVals.end() ; ++valIter) {
      file << (*valIter) << endl;
    }
  file << "        </DataArray>\n";
  file << "      </CellData>\n";
	file << "    </Piece>\n";
	file << "  </UnstructuredGrid>\n";
	file << "</VTKFile>\n";

}



void DofEdgeVectorPD::interpol(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIterP = edgeVals.begin();
  vector<double>::iterator valIterD = edgeDualVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIterP != edgeVals.end() || valIterD != edgeDualVals.end() ; 
                ++edgeIter, ++valIterP, ++valIterD) {
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
    WorldVector<double> ce = 0.5 * (p + q); // circumcenter of the primal edge
    WorldVector<double> e = q - p; // primal edge = q - p
    // set primal val
    (*valIterP) = (*alpha)(ce, e);
    
    double lenP = std::sqrt(e*e); // length of primal edge 
    p = edgeIter->infoRight->getCircumcenter();
    q = edgeIter->infoLeft->getCircumcenter();
    e = q - p; // dual edge
    // double lenD = sqrt(e*e); // length of dual edge // see below
    double lenD =  edgeIter->infoRight->getDualEdgeLen(edgeIter->dofEdge)
                  +edgeIter->infoLeft->getDualEdgeLen(edgeIter->dofEdge); // length of dual edge
    // use the identity <*alpha, primal edge> = - (lenP / lenD) <alpha, dual edge>
    (*valIterD) = - lenP * (*alpha)(ce, e) / lenD;
  }
}

void DofEdgeVectorPD::interpol(AbstractFunction<WorldVector<double>, WorldVector<double> > *vec) {
  DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
  edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIterP = edgeVals.begin();
  vector<double>::iterator valIterD = edgeDualVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end() || valIterP != edgeVals.end() || valIterD != edgeDualVals.end() ; 
                ++edgeIter, ++valIterP, ++valIterD) {
    WorldVector<double> p = coords[edgeIter->dofEdge.first];
    WorldVector<double> q = coords[edgeIter->dofEdge.second];
    WorldVector<double> ce = 0.5 * (p + q); // circumcenter of the primal edge
    WorldVector<double> vecval = (*vec)(ce);
    WorldVector<double> e = q - p; // primal edge = q - p
    // set primal val
    (*valIterP) = vecval * e;
    
    double lenP = std::sqrt(e*e); // length of primal edge 
    p = edgeIter->infoRight->getCircumcenter();
    q = edgeIter->infoLeft->getCircumcenter();
    e = q - p; // dual edge
    //double lenD = sqrt(e*e); // length of dual edge //WRONG! (Underestimated for Curvature not 0, but -> 0 for h->0)
    double lenD =  edgeIter->infoRight->getDualEdgeLen(edgeIter->dofEdge)
                  +edgeIter->infoLeft->getDualEdgeLen(edgeIter->dofEdge); // length of dual edge
    // use the identity <*alpha, primal edge> = - (lenP / lenD) <alpha, dual edge>
    (*valIterD) = - lenP * (vecval * e) / lenD;
    //(*valIterD) = - lenP * ( ((*vec)(p))*(ce - p) + ((*vec)(q))*(q - ce) ) / lenD;
  }
}

void DofEdgeVectorPD::normalize(double eps) {
  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIterP = edgeVals.begin();
  vector<double>::iterator valIterD = edgeDualVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter, ++valIterP, ++valIterD) {
    double lenP = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    double norm = std::sqrt((*valIterP) * (*valIterP) + (*valIterD) * (*valIterD)) / lenP;
    double normInv = (norm < eps) ? 0.0 : (1.0/norm);
    (*valIterP) *= normInv;
    (*valIterD) *= normInv;
  }
}

DofEdgeVector DofEdgeVectorPD::getNormOnEdges() const{
  DofEdgeVector norms(edgeMesh, name + "Magnitude");

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::const_iterator valIterP = edgeVals.begin();
  vector<double>::const_iterator valIterD = edgeDualVals.begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter, ++valIterP, ++valIterD) {
    double lenP = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    norms[edgeIter->edgeDof] = std::sqrt((*valIterP) * (*valIterP) + (*valIterD) * (*valIterD)) / lenP;
  }

  return norms;
}

double DofEdgeVectorPD::L2Norm2(const DofEdgeVector &primal, const DofEdgeVector &dual){
  DofEdgeVector norms(primal.getEdgeMesh(), "We've got the PMA!");

  vector<EdgeElement>::const_iterator edgeIter = primal.getEdgeMesh()->getEdges()->begin();
  vector<double>::const_iterator valIterP = primal.getEdgeVector()->begin();
  vector<double>::const_iterator valIterD = dual.getEdgeVector()->begin();
  for (; edgeIter != primal.getEdgeMesh()->getEdges()->end(); ++edgeIter, ++valIterP, ++valIterD) {
    double lenP = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    norms[edgeIter->edgeDof] = ((*valIterP) * (*valIterP) + (*valIterD) * (*valIterD)) / (lenP * lenP);
  }
  
  return norms.surfaceIntegration();
}

double DofEdgeVectorPD::L2Inner(const DofEdgeVector &primalA, const DofEdgeVector &dualA, const DofEdgeVector &primalB, const DofEdgeVector &dualB){
  DofEdgeVector inner(primalA.getEdgeMesh(), "We've got the PMA!");

  vector<EdgeElement>::const_iterator edgeIter = primalA.getEdgeMesh()->getEdges()->begin();
  vector<double>::const_iterator valIterPA = primalA.getEdgeVector()->begin();
  vector<double>::const_iterator valIterDA = dualA.getEdgeVector()->begin();
  vector<double>::const_iterator valIterPB = primalB.getEdgeVector()->begin();
  vector<double>::const_iterator valIterDB = dualB.getEdgeVector()->begin();
  for (; edgeIter != primalA.getEdgeMesh()->getEdges()->end(); ++edgeIter, ++valIterPA, ++valIterDA, ++valIterPB, ++valIterDB) {
    double lenP = edgeIter->infoLeft->getEdgeLen(edgeIter->dofEdge);
    inner[edgeIter->edgeDof] = ((*valIterPA) * (*valIterPB) + (*valIterDA) * (*valIterDB)) / (lenP * lenP);
  }
  
  return inner.surfaceIntegration();
}

vector<WorldVector<double> > DofEdgeVectorPD::getSharpVecOnEdges() {
  vector<WorldVector<double> > sharp(edgeMesh->getNumberOfEdges());

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  vector<double>::iterator valIterP = edgeVals.begin();
  vector<double>::iterator valIterD = edgeDualVals.begin();
  vector<WorldVector<double> >::iterator sharpIter = sharp.begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter, ++valIterP, ++valIterD, ++sharpIter) {
    DofEdge dofe = edgeIter->dofEdge;
    double lenP = edgeIter->infoLeft->getEdgeLen(dofe);
    double lenD = edgeIter->infoLeft->getDualEdgeLen(dofe) + edgeIter->infoRight->getDualEdgeLen(dofe);
    WorldVector<double> eP = edgeIter->infoLeft->getEdge(dofe);
    WorldVector<double> eD = edgeIter->infoLeft->getCircumcenter() - edgeIter->infoRight->getCircumcenter();
    (*sharpIter) = (1.0 / lenP) * ( ((*valIterP) / lenP) * eP - ((*valIterD) / lenD) * eD);
  }

  return sharp;
}

DofVertexVector DofEdgeVectorPD::interiorProdOnVertices(const DofEdgeVectorPD &pdvec) const {
  DofVertexVector iprod(edgeMesh, name + " interprod " + pdvec.getName());
  iprod.set(0.0);
  
  DOFVector<bool> visited(edgeMesh->getFeSpace(),"visited");
  visited.set(false);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    //first verex
    DegreeOfFreedom vdof = edgeIter->dofEdge.first;
    if (!visited[vdof]) {
      double vvol4 = 0.0;
      for (EdgeElement::EdgeRingIterator ringIter(&(*edgeIter), FIRSTVERTEX); !ringIter.isEnd(); ++ringIter) {
        double lenP = ringIter->infoLeft->getEdgeLen(ringIter->dofEdge);
        double lenD = ringIter->infoLeft->getDualEdgeLen(ringIter->dofEdge)
                     +ringIter->infoRight->getDualEdgeLen(ringIter->dofEdge);
        vvol4 += lenP * lenD;
        iprod[vdof] += (lenD / lenP) * ( getPrimalVal(*ringIter) * pdvec.getPrimalVal(*ringIter) 
                                      + getDualVal(*ringIter) * pdvec.getDualVal(*ringIter)    );
      }
      iprod[vdof] /=  vvol4;
      visited[vdof] = true;
    }
    //second vertex
    vdof = edgeIter->dofEdge.second;
    if (!visited[vdof]) {
      double vvol4 = 0.0;
      for (EdgeElement::EdgeRingIterator ringIter(&(*edgeIter), SECONDVERTEX); !ringIter.isEnd(); ++ringIter) {
        double lenP = ringIter->infoLeft->getEdgeLen(ringIter->dofEdge);
        double lenD = ringIter->infoLeft->getDualEdgeLen(ringIter->dofEdge)
                     +ringIter->infoRight->getDualEdgeLen(ringIter->dofEdge);
        vvol4 += lenP * lenD;
        iprod[vdof] += (lenD / lenP) *  ( getPrimalVal(*ringIter) * pdvec.getPrimalVal(*ringIter) 
                                      + getDualVal(*ringIter) * pdvec.getDualVal(*ringIter)    );
      }
      iprod[vdof] /=  vvol4;
      visited[vdof] = true;
    }
  }

  return iprod;
}

DofVertexVector DofEdgeVectorPD::interiorProdOnVertices(const DofEdgeVector &primalA, const DofEdgeVector &dualA, const DofEdgeVector &primalB, const DofEdgeVector &dualB) {
  const EdgeMesh *edgeMesh = primalA.getEdgeMesh();
  DofVertexVector iprod(edgeMesh, primalA.getName() + "_interprod_" + primalB.getName());
  iprod.set(0.0);
  
  DOFVector<bool> visited(edgeMesh->getFeSpace(),"visited");
  visited.set(false);

  vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
  for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
    //first verex
    DegreeOfFreedom vdof = edgeIter->dofEdge.first;
    if (!visited[vdof]) {
      double vvol4 = 0.0;
      for (EdgeElement::EdgeRingIterator ringIter(&(*edgeIter), FIRSTVERTEX); !ringIter.isEnd(); ++ringIter) {
        double lenP = ringIter->infoLeft->getEdgeLen(ringIter->dofEdge);
        double lenD = ringIter->infoLeft->getDualEdgeLen(ringIter->dofEdge)
                     +ringIter->infoRight->getDualEdgeLen(ringIter->dofEdge);
        vvol4 += lenP * lenD;
        iprod[vdof] += (lenD / lenP) * ( primalA[*ringIter] * primalB[*ringIter]
                                      + dualA[*ringIter] * dualB[*ringIter]    );
      }
      iprod[vdof] /=  vvol4;
      visited[vdof] = true;
    }
    //second vertex
    vdof = edgeIter->dofEdge.second;
    if (!visited[vdof]) {
      double vvol4 = 0.0;
      for (EdgeElement::EdgeRingIterator ringIter(&(*edgeIter), SECONDVERTEX); !ringIter.isEnd(); ++ringIter) {
        double lenP = ringIter->infoLeft->getEdgeLen(ringIter->dofEdge);
        double lenD = ringIter->infoLeft->getDualEdgeLen(ringIter->dofEdge)
                     +ringIter->infoRight->getDualEdgeLen(ringIter->dofEdge);
        vvol4 += lenP * lenD;
        iprod[vdof] += (lenD / lenP) *  ( primalA[*ringIter] * primalB[*ringIter]
                                      + dualA[*ringIter] * dualB[*ringIter]    );
      }
      iprod[vdof] /=  vvol4;
      visited[vdof] = true;
    }
  }

  return iprod;
}


void DofEdgeVectorPD::writeSharpOnEdgesFile(string name) {
    boost::iostreams::filtering_ostream file;
    file.push(boost::iostreams::file_descriptor_sink(name, std::ios::trunc));
  file << "<?xml version=\"1.0\"?>\n";
	file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	file << "  <UnstructuredGrid>\n";
	file << "    <Piece NumberOfPoints=\"" << edgeMesh->getFeSpace()->getMesh()->getNumberOfVertices()
	     << "\" NumberOfCells=\"" << edgeVals.size() << "\">\n";
	file << "      <Points>\n";
	file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    DOFVector< WorldVector< double > > coords(edgeMesh->getFeSpace(), "coords");
    edgeMesh->getFeSpace()->getMesh()->getDofIndexCoords(coords);
    DOFAdmin *admin = edgeMesh->getFeSpace()->getAdmin();
    for (DegreeOfFreedom dof = 0; dof < coords.getSize() ; dof++) {
      if (!admin->isDofFree(dof)) {
        file << coords[dof] << endl;
      }
    }
  file << "        </DataArray>\n";
	file << "      </Points>\n";
	file << "      <Cells>\n";
	file << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";
    for (int i = 1; i <= edgeVals.size(); ++i) {
      file << (i *  2) << endl;
    }
  file << "        </DataArray>\n";
	file << "        <DataArray type=\"UInt8\" Name=\"types\">\n";
    for (int i = 1; i <= edgeVals.size(); ++i) {
      file << 3 << endl;
    }
  file << "        </DataArray>\n";
	file << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";
    vector<EdgeElement>::const_iterator edgeIter = edgeMesh->getEdges()->begin();
    for (; edgeIter != edgeMesh->getEdges()->end(); ++edgeIter) {
      file << edgeIter->dofEdge.first << " " << edgeIter->dofEdge.second << endl;
    }
  file << "        </DataArray>\n";    
	file << "      </Cells>\n";
	file << "      <CellData>\n";
  file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\""
	     << this->name
	     << "\" format=\"ascii\">\n";
    vector<WorldVector<double> > sharp = getSharpVecOnEdges();
    vector<WorldVector<double> >::iterator sharpIter = sharp.begin();
    for (; sharpIter != sharp.end() ; ++sharpIter) {
      file << (*sharpIter) << endl;
    }
  file << "        </DataArray>\n";
  file << "      </CellData>\n";
	file << "    </Piece>\n";
	file << "  </UnstructuredGrid>\n";
	file << "</VTKFile>\n";

}
