#include "AMDiS.h"

namespace AMDiS {

DOFVector<int> getConnections(const FiniteElemSpace *feSpace);

DOFVector<double> get1RingVols(const FiniteElemSpace *feSpace);

DOFVector<double> getDualVols(const FiniteElemSpace *feSpace);

// 1-Ring approximation, don't need wellcentered
DOFVector<double> getVoronoiRadii(const FiniteElemSpace *feSpace);

// need wellcentered
DOFVector<double> getVoronoiRadiiDualApprox(const FiniteElemSpace *feSpace);

// weighted, not nomalized
DOFVector<WorldVector<double> > getNormals(const FiniteElemSpace *feSpace);

DOFVector<WorldVector<double> > getConnectionForces(const FiniteElemSpace *feSpace, bool constantRadii = false, double k = 0.75);

double getMaxMagnitude(DOFVector<WorldVector<double> > F);

  
}
