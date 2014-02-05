#include "AMDiS.h"

namespace AMDiS {

DOFVector<int> getConnections(const FiniteElemSpace *feSpace);

DOFVector<double> get1RingVols(const FiniteElemSpace *feSpace);

DOFVector<double> getDualVols(const FiniteElemSpace *feSpace);

// 1-Ring approximation, don't need wellcentered
DOFVector<double> getVoronoiRadii(const FiniteElemSpace *feSpace);

// need wellcentered
DOFVector<double> getVoronoiRadiiDualApprox(const FiniteElemSpace *feSpace);

// k...scale parameter
DOFVector<WorldVector<double> > getConnectionForces(const FiniteElemSpace *feSpace, bool constantRadii = false, double k = 1.0);
  
}
