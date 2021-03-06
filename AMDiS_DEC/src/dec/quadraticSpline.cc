#include "quadraticSpline.h"
#include "DOFVHelper.h"
#include "MeshHelper.h"
#include "WorldVectorHelper.h"

QuadraticSpline::QuadraticSpline(const FiniteElemSpace *finiteElemSpace) : 
    feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace)),
    coeffs(feSpace, "coeffs")
{
  cout << "*** calculate quardratic spline surface approximation ***" << endl;
  coeffs = mtl::dense_vector<double>(10);

  DOFVector<mtl::dense2D<double> > terms(feSpace, "qwertz");
  terms = mtl::dense2D<double>(0, 10);
  
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
      WorldVector<double> x = el->getCoord(i);
      WorldVector<double> next_x = el->getCoord((i+1)%3) - x;
      WorldVector<double> opp_x = el->getOppCoord(i) - x;
      // TODO: expensive?
      if (visited[dof]) {
        locMat->change_dim(r + 2, 10, true);
      } else {
        locMat->change_dim(r + 3, 10, true);
      }
      for (int k = 0; k < 3; k++) {
        int k1 = (k+1)%3;
        int k2 = (k+2)%3;
      // only next vertex by orientation to prevent double-entry
        (*locMat)[r][k] = 0.5 * next_x[k] * next_x[k];
        (*locMat)[r][k+3] = next_x[k1] * next_x[k2];
        (*locMat)[r][k+6]= next_x[k];
        (*locMat)[r][9]= 1.0;
        // opposite vertex
        (*locMat)[r+1][k] = 0.5 * opp_x[k] * opp_x[k];
        (*locMat)[r+1][k+3] = opp_x[k1] * opp_x[k2];
        (*locMat)[r+1][k+6]= opp_x[k];
        (*locMat)[r+1][9]= 1.0;
      }
      if (!visited[dof]) {
        for (int k = 0; k < 3; k++){
        // vertex self
          (*locMat)[r+2][k] = 0.5 * x[k] * x[k];
          (*locMat)[r+2][k+3] = x[(k+1)%3] * x[(k+2)%3];
          (*locMat)[r+2][k+6]= x[k];
          (*locMat)[r+2][9]= 1.0;
        }
        visited[dof] = true;
      }
    }
  }
  cout << "*** build terms done ***" << endl;

  DOFVector<mtl::dense2D<double> >::Iterator termsIter(&terms, USED_DOFS);
  DOFVector<mtl::dense_vector<double> >::Iterator coeffsIter(&coeffs, USED_DOFS);
  
  for (termsIter.reset(), coeffsIter.reset(); !coeffsIter.end(); ++termsIter, ++coeffsIter) {
    mtl::dense2D<double> sysMat(10,10); 
    sysMat = trans(*termsIter) * (*termsIter);
    mtl::dense_vector<double> rhs(10);
    rhs = 0.0;
    //rhs = trans(*termsIter) * mtl::dense_vector<double>(num_rows(*termsIter), 1.0);
    
    ////precond
    //mtl::dense2D<double> L(9,9);
    //L = 0.0;
    //for (int i = 0 ; i < 9; i++) {
    //  L[i][i] = 1./sysMat[i][i];
    //}
    //mtl::dense2D<double> LA(9,9);
    //LA = L * sysMat;
    //mtl::dense_vector<double> Lb(9);
    //Lb = L * rhs;

    //*coeffsIter = inv(LA) * Lb;
    *coeffsIter = inv(sysMat) * rhs;
    mtl::dense_vector<double> eig(10) ;
    //eig = qr_algo(LA,100);
    eig = qr_algo(sysMat,100);
    cout << (eig[0]/eig[9]) << endl;
    //cout << (*termsIter) << endl;
    //cout << LA << endl;

    // precond
    //itl::pc::diagonal<mtl::dense2D<double> > L(sysMat), R(sysMat);
    //*coeffsIter = 0.0;
    //// Termination criterion: r < 1e-6 * b or N iterations
    //itl::noisy_iteration<double>       iter(rhs, 1000, 1.e-6);
    //iter.set_quite(true);
    //iter.suppress_resume(true);
    //itl::tfqmr(sysMat, *coeffsIter, rhs, L, R, iter);

  }
  cout << "*** solve systems done ***" << endl;

}

DOFVector<WorldVector<double> > QuadraticSpline::getNormals() 
{
  DOFVector<WorldVector<double> > nu(feSpace, "normals");
  //for orientation -- improve
  DOFVector<WorldVector<double> > normals = AMDiS::getNormals(feSpace);

  const BasisFunction *basFcts = feSpace->getBasisFcts();
  int numBasFcts = basFcts->getNumber();
  std::vector<DegreeOfFreedom> localIndices(numBasFcts);
  DOFAdmin *admin = feSpace->getAdmin();
  DegreeOfFreedom dof;

  std::map<DegreeOfFreedom, bool> visited;
  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS); el; el = stack.traverseNext(el)) 
  {
    basFcts->getLocalIndices(el->getElement(), admin, localIndices);
    for (int i = 0; i < 3; i++) {
      dof = localIndices[i];
      if (!visited[dof]) {
    //cout << coeffs[dof] << endl;
        WorldVector<double> x = el->getCoord(i);
        for (int k = 0; k < 3; k++){
          int k1 = (k+1)%3;
          int k2 = (k+2)%3;
          nu[dof][k] = coeffs[dof][k] * x[k] + coeffs[dof][6-k-k1] * x[k1] + coeffs[dof][6-k-k2] * x[k2] + coeffs[dof][6+k];
        }
  //for orientation -- improve
        if (dot(normals[dof], nu[dof]) <= 0) nu[dof] *= -1.0;
        visited[dof] = true;
      }
    }
  }

  return normalize(nu);
}
