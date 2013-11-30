#include "AMDiS.h"

namespace AMDiS {

class ElVolumesInfo2d {
  
  public:
    ElVolumesInfo2d(ElInfo *el);

  protected:
    WorldVector<double> dualVertexVol;
    WorldVector<double> oppEdgeLen;
    WorldVector<double> dualOppEdgeLen;
};
}
