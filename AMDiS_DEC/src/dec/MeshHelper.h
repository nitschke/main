#include "AMDiS.h"

namespace AMDiS {

DOFVector<int> getConnections(const FiniteElemSpace *feSpace);

DOFVector<double> get1RingVols(const FiniteElemSpace *feSpace);

DOFVector<double> getDualVols(const FiniteElemSpace *feSpace);

DOFVector<double> getVoronoiRadii(const FiniteElemSpace *feSpace);
  
}
