#include "quadratic.h"
#include "DOFVHelper.h"
#include "phiProjection.h"

Quadratic::Quadratic(const FiniteElemSpace *finiteElemSpace) : 
  coeffs(9),
  feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace))
{
  cout << "*** calculate quardratic surface approximation ***" << endl;
  mtl::dense2D<double> mat(feSpace->getMesh()->getNumberOfVertices(), 9);
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

  mtl::dense2D<double> sysMat(9,9); 
  sysMat = trans(mat) * mat;
  mtl::dense_vector<double> rhs(9); 
  rhs = trans(mat) * mtl::dense_vector<double>(feSpace->getMesh()->getNumberOfVertices(), 1.0);
  
  coeffs = inv(sysMat) * rhs;
  cout << "*** done ***" << endl;
}

DOFVector<WorldVector<double> > Quadratic::getNormals() 
{
  DOFVector<WorldVector<double> > nu(feSpace, "normals");

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
          int k1 = (k+1)%3;
          int k2 = (k+2)%3;
          nu[dof][k] = coeffs[k] * x[k] + coeffs[6-k-k1] * x[k1] + coeffs[6-k-k2] * x[k2] + coeffs[6+k];
        }
        visited[dof] = true;
      }
    }
  }

  return normalize(nu);
}

void Quadratic::writeVTK(string name)
{
  new PhiProject(42, VOLUME_PROJECTION, new Phi(coeffs), new GradPhi(coeffs), 1.0e-8);

  WorldVector<DOFVector<double> * > coords;
  for (int i = 0; i < 3; i++) 
  {
      coords[i] = new DOFVector<double>(feSpace, "paraCoordComponent");
  }
  const BasisFunction *basFcts = feSpace->getBasisFcts();
  int numBasFcts = basFcts->getNumber();
  std::vector<DegreeOfFreedom> localIndices(numBasFcts);
  DOFAdmin *admin = feSpace->getAdmin();
  DegreeOfFreedom dof;

  std::map<DegreeOfFreedom, bool> visited;
  TraverseStack stack;
  for (ElInfo *el = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS); el; el = stack.traverseNext(el)) {
    basFcts->getLocalIndices(el->getElement(), admin, localIndices);
    for (int i = 0; i < 3; i++) {
      dof = localIndices[i];
      if (!visited[dof]) {
        WorldVector<double> x =  el->getCoord(i);
        Projection::getProjection(42)->project(x);
        for (int k = 0; k < 3; k++){
          (*(coords[k]))[dof] = x[k];
        }
        visited[dof] = true;
      }
    }
  }
  Parametric *parametric = feSpace->getMesh()->getParametric(); 
  if (parametric) delete parametric;
  parametric = new ParametricFirstOrder(&coords);
  feSpace->getMesh()->setParametric(parametric);

  DOFVector<double> dummy(feSpace,"dummy");
  VtkVectorWriter::writeFile(dummy, name);
}
