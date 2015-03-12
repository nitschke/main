#ifndef ELVOLUMESINFO2D_H
#define ELVOLUMESINFO2D_H

#include "AMDiS.h"

namespace AMDiS {

class ElVolumesInfo2d {
  
  public:
    ElVolumesInfo2d(const ElInfo *el);
    
    const ElInfo *getElInfo() const {return elInfo;}

    double getVol() {return 0.5 * elInfo->getDet();}

    int getLocal(DegreeOfFreedom dof) {
      std::map<DegreeOfFreedom, int>::iterator it = dofToLocal.find(dof);
      TEST_EXIT(it != dofToLocal.end())("ElVolumesInfo2d::getLocal: Element dont contain this global DOF!");
      return it->second;
    }

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

    std::map<DegreeOfFreedom, int> dofToLocal;

    const ElInfo *elInfo;
};
}

#endif
