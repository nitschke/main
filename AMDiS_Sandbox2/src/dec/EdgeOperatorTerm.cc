#include "EdgeOperatorTerm.h"

edgeRowValMapper Discrete1FormAtEdge::evalRow(EdgeElement *eel, double factor) {
  edgeRowValMapper rowMapper;
  DegreeOfFreedom edof = eel->edgeDof;
  rowMapper[edof] = fac * factor * evec[edof];
  return rowMapper;
}
