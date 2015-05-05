#include "AMDiS.h"
#include "torusProjection.h"

namespace AMDiS {

  void TorusProject::project(WorldVector<double> &x) 
  {
    WorldVector<double> uTilde= x;
    uTilde[1]= 0.0;;
    WorldVector<double> u= uTilde*(R/norm(uTilde));
    WorldVector<double> xTilde= x - u;
    WorldVector<double> yTilde= xTilde*(r/norm(xTilde));
    x= u + yTilde;
  }

}
