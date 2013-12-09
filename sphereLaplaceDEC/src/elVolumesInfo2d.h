#include "AMDiS.h"

namespace AMDiS {

class ElVolumesInfo2d {
  
  public:
    ElVolumesInfo2d(ElInfo *el);
    
    double getDualVertexVol(int i) {return dualVertexVol[i]}
    double getOppEdgeLen(int i) {return oppEdgeLen[i]}
    double getDualOppEdgeLen(int i) {return dualOppEdgeLen[i]}

  protected:
    WorldVector<double> dualVertexVol;
    WorldVector<double> oppEdgeLen;
    WorldVector<double> dualOppEdgeLen;
};
}
