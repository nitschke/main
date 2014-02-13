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

  DOFVector<double> Vols = get1RingVols(feSpace);
  
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
      Radii[dof] += 0.5 * el->getDet() * elR[i] / Vols[dof];
    }
  }

  return Radii;
}


DOFVector<double> getVoronoiRadiiDualApprox(const FiniteElemSpace *feSpace) {
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

DOFVector<WorldVector<double> > getNormals(const FiniteElemSpace *feSpace) {
  
  DOFVector<WorldVector<double> > normals(feSpace, "normals");
  WorldVector<double> zero(DEFAULT_VALUE, 0.0);
  normals = zero;

  DOFVector<double> vol1Ring = get1RingVols(feSpace);
  
  const BasisFunction *basFcts = feSpace->getBasisFcts();
    int numBasFcts = basFcts->getNumber();
    std::vector<DegreeOfFreedom> localIndices(numBasFcts);
    DOFAdmin *admin = feSpace->getAdmin();
    DegreeOfFreedom dof;

  WorldVector<double> normal;

  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET); el; el = stack.traverseNext(el)) {
    basFcts->getLocalIndices(el->getElement(), admin, localIndices);
    el->getElementNormal(normal);
    for (int i = 0; i < 3; i++) {
      dof = localIndices[i];
      normals[dof] += 0.5 * (el->getDet()/vol1Ring[dof]) * normal;
    }
  }

  return normals;
}

DOFVector<WorldVector<double> > getConnectionForces(const FiniteElemSpace *feSpace, bool constantRadii, double k) {
  DOFVector<WorldVector<double> > F(feSpace, "ConnectionForces");
  WorldVector<double> zero(DEFAULT_VALUE, 0.0);
  F = zero;

  double kInit = -1.0;
  Parameters::get("edgeForces->k", kInit);
  if (kInit != -1.0) k = kInit;

  DOFVector<double> Radii = getVoronoiRadii(feSpace);
  double lRef = 2.0 * Radii.average();
  //if (constantRadii) Radii = 0.5 * Radii.average();
  //Radii = 0.0;

  DOFVector<WorldVector<double> > normals = getNormals(feSpace);

  const BasisFunction *basFcts = feSpace->getBasisFcts();
    int numBasFcts = basFcts->getNumber();
    std::vector<DegreeOfFreedom> localIndices(numBasFcts);
    DOFAdmin *admin = feSpace->getAdmin();
    DegreeOfFreedom dofi;
    DegreeOfFreedom dofj;

  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET); el; el = stack.traverseNext(el)) {
    basFcts->getLocalIndices(el->getElement(), admin, localIndices);
    //ElVolumesInfo2d volInfo(el);
    for (int i = 0; i < 3; i++) {
      int j = (i+1)%3;
      int m = (i+1)%3;
      dofi = localIndices[i];
      dofj = localIndices[j];

      //WorldVector<double> ccjk = el->getCoord(j) + 0.5 * (el->getCoord(k) - el->getCoord(j));
      //WorldVector<double> ccccjk = ccjk - volInfo.getCircumcenter();
      //double len2 = sqrt(dot(ccccjk, ccccjk));
      WorldVector<double> xDelta = el->getCoord(j) - el->getCoord(i);
      WorldVector<double> xDeltaK = el->getCoord(m) - el->getCoord(i);
      double len = sqrt(dot(xDelta, xDelta));
      double lenK = sqrt(dot(xDeltaK, xDeltaK));
      double cosAngle = dot(xDelta, xDeltaK) / len / lenK;
      double deltaLen = len / lRef - k;
      double deltaA = cosAngle - 0.3;
      //linear law
      double c = 1./2.;
      double fe = c * deltaLen + (1.0-c) * deltaA;
      WorldVector<double> Fe = (fe/len) * xDelta;
      F[dofi] += Fe - (dot(Fe, normals[dofi]) / dot(normals[dofi],normals[dofi])) * normals[dofi];
      F[dofj] -= Fe - (dot(Fe, normals[dofj]) / dot(normals[dofj],normals[dofj])) * normals[dofj];
    }
  }

  return F;
}

double getMaxMagnitude(DOFVector<WorldVector<double> > F) {
  DOFIterator<WorldVector<double> > fIter(&F, USED_DOFS);
  double max = 0;
  for (fIter.reset(); !fIter.end(); fIter++) {
    double mag = sqrt(dot(*fIter,*fIter));
    if (mag > max) max = mag;
  }
  return max;
}

DOFVector<double> getAverage(DOFVector<double> f) {
  const FiniteElemSpace *feSpace = f.getFeSpace();
  DOFVector<double> rval(feSpace, f.getName() + "Average");
  rval = 0.0;

  DOFVector<double> vols = getDualVols(feSpace);

  const BasisFunction *basFcts = feSpace->getBasisFcts();
  int numBasFcts = basFcts->getNumber();
  std::vector<DegreeOfFreedom> localIndices(numBasFcts);
  DOFAdmin *admin = feSpace->getAdmin();
  DegreeOfFreedom dofi;
  DegreeOfFreedom dofj;

  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET); el; el = stack.traverseNext(el)) {
    basFcts->getLocalIndices(el->getElement(), admin, localIndices);
    
    for (int i = 0; i < 3; i++) {
      dofi = localIndices[i];

      ElVolumesInfo2d elVolInfo(el);
      double c = 0.5 * elVolInfo.getDualVertexVol(i) / vols[dofi];
      rval[dofi] += c * f[dofi];
      
      c *= 0.5;
      for (int j = (i+1)%3; i != j; j = (j+1)%3) {
        dofj = localIndices[j];
        rval[dofi] += c * f[dofj]; 
      }
    }
  }

  return rval;
}



void MeshInfoCSVWriter::appendData(const FiniteElemSpace *feSpace) {
  DOFVector<double> one(feSpace,"jhsdfakf");
  one = 1.0;
  double volM = one.Int();
  //one.~DOFVector<double>();

  double avDia = 0.0;
  double maxDia = 0.0;
  double avArea = 0.0;
  double maxArea = 0.0;
  double avAngle = 0.0;
  double maxAngle = 0.0;
  double avRat = 0.0;
  double maxRat = 0.0; 


  TraverseStack stack;
  ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET);
  for (; el; el = stack.traverseNext(el)) {
    double vol = 0.5 * el->getDet();
    double c = vol / volM;

    ElVolumesInfo2d volInfo(el);
    double dia = volInfo.getDiameter();
    if (dia > maxDia) maxDia = dia;
    avDia += c * dia;

    if (vol > maxArea) maxArea = vol;
    avArea += c * vol;

    Vector<WorldVector<double> > edges(3);
    for (int i = 0; i < 3; i++) {
      edges[i] = el->getCoord((i+1)%3) - el->getCoord(i);
    }
    double amax = 0.0;
    double amin = 4.0;
    for (int i = 0; i < 3; i++) {
      int ii = (i+1)%3;
      double a = acos(dot(edges[i], (-1.0)*edges[ii]) / wvnorm(edges[i]) / wvnorm(edges[ii]));
      if (a > amax) amax = a;
      if (a < amin) amin = a;
    }
    if (amax > maxAngle) maxAngle = amax;
    avAngle += c * amax;
    double arat = amax / amin;
    if (arat > maxRat) maxRat = arat;
    avRat += c * arat;
  }
  //stack.~TraverseStack();

  out << n++;
  out << "," << avDia;
  out << "," << maxDia;
  out << "," << avArea;
  out << "," << maxArea;
  out << "," << avAngle * 180.0 / M_PI;
  out << "," << maxAngle * 180.0 / M_PI ;
  out << "," << avRat;
  out << "," << maxRat << endl;

}

}
