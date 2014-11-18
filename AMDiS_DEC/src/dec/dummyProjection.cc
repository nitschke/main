#include "AMDiS.h"
#include "dummyProjection.h"
#include "WorldVectorHelper.h"

namespace AMDiS {

  void DummyProject::project(WorldVector<double> &x) {
  }

  void DummyProject::project(const WorldVector<double> &x, WorldVector<double> &v) {
  }

}
