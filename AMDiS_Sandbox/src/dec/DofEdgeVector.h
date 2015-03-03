#include "AMDiS.h"
#include "EdgeMesh.h"

using namespace AMDiS;
using namespace std;

class DofEdgeVector {

  DofEdgeVector(const EdgeMesh *edgeMesh_, std:string name_) : 
      edgeMesh(edgeMesh_), name(name_), edgeVals(edgeMesh->getNumberOfEdges()) {}

  // Gauss-Legendre-Intigration (n=4) of I_edge[alpha]
  void interpolGL4(BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> > *alpha,
            AbstractFunction<WorldVector<double>, WorldVector<double> > *proj,
            AbstractFunction<WorldMatrix<double>, WorldVector<double> > *jproj);

  
 
private:
  const EdgeMesh *edgeMesh;

  vector<double> edgeVals;

  std::string name;

};
