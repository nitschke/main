#include "AMDiS.h"

namespace AMDiS {

class ElVolumesInfo2d {
  
  public:
    ElVolumesInfo2d(const ElInfo *el);
    
    double getDualVertexVol(int i) {return dualVertexVol[i];}
    double getOppEdgeLen(int i) {return oppEdgeLen[i];}
    double getDualOppEdgeLen(int i) {return dualOppEdgeLen[i];}
    WorldVector<double> getCircumcenter() {return cc;}

    double getDiameter();

    // Angle on vertex
    double getAngle(int i);

    // Sin(Angle) on vertex
    double getSin(int i);

  protected:
    WorldVector<double> dualVertexVol;
    WorldVector<double> oppEdgeLen;
    WorldVector<double> dualOppEdgeLen;
    WorldVector<double> cc;

    const ElInfo *elInfo;
};
}
