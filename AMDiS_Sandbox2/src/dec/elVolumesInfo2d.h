#ifndef ELVOLUMESINFO2D_H
#define ELVOLUMESINFO2D_H

#include "AMDiS.h"

using namespace AMDiS;
using namespace std;


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

    int getOppVertexLocal(DofEdge dofEdge) {
      return 3 - getLocal(dofEdge.first) - getLocal(dofEdge.second);
    }

    double getEdgeLen(DofEdge dofEdge) {
      return getOppEdgeLen(getOppVertexLocal(dofEdge));
    }

    double getDualEdgeLen(DofEdge dofEdge) {
      return getDualOppEdgeLen(getOppVertexLocal(dofEdge));
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
    
inline ostream &operator <<(ostream &out, const ElVolumesInfo2d &volInfo) {
  Element *el = volInfo.getElInfo()->getElement();
  for (int i = 0; i < 3 ; ++i) {
    out << el->getDof(i,0) << " "; 
  }
  return out;
}





#endif
