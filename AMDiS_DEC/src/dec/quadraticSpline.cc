#include "quadraticSpline.h"

QuadraticSpline::QuadraticSpline(const FiniteElemSpace *finiteElemSpace) : 
    feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace)),
    coeffs(feSpace, "coeffs")
{
  cout << "*** calculate quardratic spline surface approximation ***" << endl;
  coeffs = mtl::dense_vector<double>(9);

  DOFVector<mtl::dense2D<double> > terms(feSpace, "qwertz");
  terms = mtl::dense2D<double>(0, 9);
  
  const BasisFunction *basFcts = feSpace->getBasisFcts();
  int numBasFcts = basFcts->getNumber();
  std::vector<DegreeOfFreedom> localIndices(numBasFcts);
  DOFAdmin *admin = feSpace->getAdmin();
  DegreeOfFreedom dof;

  std::map<DegreeOfFreedom, bool> visited;
  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_OPP_COORDS); el; el = stack.traverseNext(el)) 
  {
    basFcts->getLocalIndices(el->getElement(), admin, localIndices);
    for (int i = 0; i < 3; i++) {
      dof = localIndices[i];
      mtl::dense2D<double> *locMat = &(terms[dof]);
      int r = num_rows(*locMat);
      WorldVector<double> next_x = el->getCoord((i+1)%3);
      WorldVector<double> opp_x = el->getOppCoord(i);
      // TODO: expensive?
      if (visited[dof]) {
        locMat->change_dim(r + 2, 9, true);
      } else {
        locMat->change_dim(r + 3, 9, true);
      }
      for (int k = 0; k < 3; k++) {
        int k1 = (k+1)%3;
        int k2 = (k+2)%3;
      // only next vertex by orientation to prevent double-entry
        (*locMat)[r][k] = 0.5 * next_x[k] * next_x[k];
        (*locMat)[r][k+3] = next_x[k1] * next_x[k2];
        (*locMat)[r][k+6]= next_x[k];
        // opposite vertex
        (*locMat)[r+1][k] = 0.5 * opp_x[k] * opp_x[k];
        (*locMat)[r+1][k+3] = opp_x[k1] * opp_x[k2];
        (*locMat)[r+1][k+6]= opp_x[k];
      }
      if (!visited[dof]) {
        WorldVector<double> x = el->getCoord(i);
        for (int k = 0; k < 3; k++){
        // vertex self
          (*locMat)[r+2][k] = 0.5 * x[k] * x[k];
          (*locMat)[r+2][k+3] = x[(k+1)%3] * x[(k+2)%3];
          (*locMat)[r+2][k+6]= x[k];
        }
        visited[dof] = true;
      }
    }
  }
  cout << "*** build terms done ***" << endl;

  DOFVector<mtl::dense2D<double> >::Iterator termsIter(&terms, USED_DOFS);
  DOFVector<mtl::dense_vector<double> >::Iterator coeffsIter(&coeffs, USED_DOFS);
  
  for (termsIter.reset(), coeffsIter.reset(); !coeffsIter.end(); ++termsIter, ++coeffsIter) {
    mtl::dense2D<double> sysMat(9,9); 
    sysMat = trans(*termsIter) * (*termsIter);
    mtl::dense_vector<double> rhs(9); 
    rhs = trans(*termsIter) * mtl::dense_vector<double>(num_rows(*termsIter), 1.0);
    //TODO: improve
    *coeffsIter = inv(sysMat) * rhs;
    cout << *coeffsIter << endl;
  }
  cout << "*** solve systems done ***" << endl;

}

DOFVector<WorldVector<double> > Quadratic::getNormals() 
{
  DOFVector<WorldVector<double> > nu(feSpace, "normals");
}
