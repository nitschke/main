#ifndef ELVOLUMESINFO2D_H
#define ELVOLUMESINFO2D_H

#include "AMDiS.h"

using namespace std;
namespace AMDiS { namespace dec {

class ElVolumesInfo2d {
  
  public:
    ElVolumesInfo2d(const ElInfo *el);
    
    const ElInfo *getElInfo() const {return elInfo;}

    double getVol() const {return 0.5 * elInfo->getDet();}

    double getDualVertexVol(int i) const {return dualVertexVol[i];}
    double getOppEdgeLen(int i) const {return oppEdgeLen[i];}
    double getDualOppEdgeLen(int i) const {return dualOppEdgeLen[i];}
    WorldVector<double> getCircumcenter() const {return cc;}

    double getDiameter();

    // Angle on vertex
    double getAngle(int i);

    // Sin(Angle) on vertex
    double getSin(int i) const;

    int getLocal(DegreeOfFreedom dof) const {
      std::map<DegreeOfFreedom, int>::const_iterator it = dofToLocal.find(dof);
      TEST_EXIT(it != dofToLocal.end())("ElVolumesInfo2d::getLocal: Element dont contain this global DOF!");
      return it->second;
    }

    WorldVector<double> getEdge(DofEdge dofEdge) {
      return elInfo->getCoord(getLocal(dofEdge.second)) - elInfo->getCoord(getLocal(dofEdge.first));
    }

    WorldVector<double> getEdgeCenter(DofEdge dofEdge) {
      return 0.5 *(elInfo->getCoord(getLocal(dofEdge.second)) + elInfo->getCoord(getLocal(dofEdge.first)));
    }

    WorldVector<double> getCoordFromGlobalIndex(DegreeOfFreedom dof) {
      return elInfo->getCoord(getLocal(dof));
    }

    int getOppVertexLocal(DofEdge dofEdge) const {
      return 3 - getLocal(dofEdge.first) - getLocal(dofEdge.second);
    }

    double getEdgeLen(DofEdge dofEdge) const {
      return getOppEdgeLen(getOppVertexLocal(dofEdge));
    }

    double getDualEdgeLen(DofEdge dofEdge) const {
      return getDualOppEdgeLen(getOppVertexLocal(dofEdge));
    }

    double getDualVertexVol_global(DegreeOfFreedom dof) const {
      return getDualVertexVol(getLocal(dof));
    }

    double getSin_global(DegreeOfFreedom dof) const {
      return getSin(getLocal(dof));
    }

    ~ElVolumesInfo2d() {
      delete elInfo;
    }



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



}}

#endif
