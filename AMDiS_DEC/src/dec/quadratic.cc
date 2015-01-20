#include "quadratic.h"

Quadratic::Quadratic(const FiniteElemSpace *finiteElemSpace) : 
  coeffs(9),
  feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace))
{
  dense2D<double> mat(feSpace->getMesh()->getNumberOfVertices(), 9);
  mat = 0.0;

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
        WorldVector<double> x = el->getCoord(i);
        for (int k = 0; k < 3; k++){
          mat[dof][k] = 0.5 * x[k] * x[k];
          mat[dof][k+3] = x[(k+1)%3] * x[(k+2)%3];
          mat[dof][k+6] = x[k];
        }
        visited[dof] = true;
      }
    }
  }

  dense2D<double> sysMat = trans(mat) * mat;
  dense_vector<double> rhs = trans(mat) * dense_vector<double>(feSpace->getMesh()->getNumberOfVertices(), 1.0);
  
  coeffs = inv(sysMat) * rhs;
  cout << coeffs << endl;
}
