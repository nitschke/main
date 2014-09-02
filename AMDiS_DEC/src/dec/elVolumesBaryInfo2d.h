#include "AMDiS.h"

namespace AMDiS {

class ElVolumesBaryInfo2d {
  
  public:
    ElVolumesBaryInfo2d(const ElInfo *el);
    
    double getDualVertexVol(int i) {return dualVertexVol[i];}
    double getOppEdgeLen(int i) {return oppEdgeLen[i];}
    double getDualOppEdgeLen(int i) {return dualOppEdgeLen[i];}
    double getOppSin(int i) {return oppSin[i];}
    WorldVector<double> getCenter() {return cc;}

    //double getDiameter();

  protected:
    WorldVector<double> dualVertexVol;
    WorldVector<double> oppEdgeLen;
    WorldVector<double> dualOppEdgeLen;
    WorldVector<double> oppSin;
    WorldVector<double> cc;

    const ElInfo *elInfo;
};
}
