#include "EdgeOperatorTerm.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"

using namespace AMDiS;
using namespace dec;


edgeRowValMapper Discrete1FormAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  DegreeOfFreedom edof = eel.edgeDof;
  rowMapper[edof] = fac * factor * (*evec)[edof];
  return rowMapper;
}

edgeRowValMapper IdentityAtEdges::evalRow(const EdgeElement &eel, double factor) {
  edgeRowValMapper rowMapper;
  rowMapper[eel.edgeDof] = fac * factor;
  return rowMapper;
}
