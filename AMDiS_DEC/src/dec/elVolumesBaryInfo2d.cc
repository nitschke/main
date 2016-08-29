#include "elVolumesBaryInfo2d.h"
#include "WorldVectorHelper.h"

namespace AMDiS {

  ElVolumesBaryInfo2d::ElVolumesBaryInfo2d(const ElInfo *el): elInfo(el) {
    dualVertexVol = WorldVector<double>();
    oppEdgeLen = WorldVector<double>();
    dualOppEdgeLen = WorldVector<double>();
    cc = WorldVector<double>();

    WorldVector<double> p0 = el->getCoord(0);
    WorldVector<double> p1 = el->getCoord(1);
    WorldVector<double> p2 = el->getCoord(2);

    WorldVector<double> e0 = p2 - p1;
    WorldVector<double> e1 = p0 - p2;
    WorldVector<double> e2 = p1 - p0;


    // barycenter
    cc = (1.0/3.0) * (p0 + p1 + p2);

    // edge circumcenters
    WorldVector<double> cc0 = p1 + 0.5*e0;
    WorldVector<double> cc1 = p2 + 0.5*e1;
    WorldVector<double> cc2 = p0 + 0.5*e2;

    // dual edges (in T)
    WorldVector<double> starE0 = cc - cc0; 
    WorldVector<double> starE1 = cc - cc1; 
    WorldVector<double> starE2 = cc - cc2; 

    // edge len of the opp vertex
    oppEdgeLen[0] = std::sqrt(dot(e0, e0));
    oppEdgeLen[1] = std::sqrt(dot(e1, e1));
    oppEdgeLen[2] = std::sqrt(dot(e2, e2));

    // dual edge len of the opp vertex
    dualOppEdgeLen[0] = std::sqrt(dot(starE0, starE0));
    dualOppEdgeLen[1] = std::sqrt(dot(starE1, starE1));
    dualOppEdgeLen[2] = std::sqrt(dot(starE2, starE2));

    // sinus of edge and dual edge
    WorldVector<double> oppCos;
    oppCos[0] = dot(e0, starE0) / (wvnorm(e0) * wvnorm(starE0));
    oppCos[1] = dot(e1, starE1) / (wvnorm(e1) * wvnorm(starE1));
    oppCos[2] = dot(e2, starE2) / (wvnorm(e2) * wvnorm(starE2));
    for (int i = 0; i < 3 ; i++) oppSin[i] = std::sqrt(1 - oppCos[i]*oppCos[i]);
    //cout << oppSin[0] << endl;
    //cout << oppSin[1] << endl;
    //cout << oppSin[2] << endl;

    // Vol of the dual vertex (voronoi cell) (in T)
    //dualVertexVol[0] = 0.25 * (dualOppEdgeLen[1]*oppEdgeLen[1] + dualOppEdgeLen[2]*oppEdgeLen[2]);
    //dualVertexVol[1] = 0.25 * (dualOppEdgeLen[0]*oppEdgeLen[0] + dualOppEdgeLen[2]*oppEdgeLen[2]);
    //dualVertexVol[2] = 0.25 * (dualOppEdgeLen[1]*oppEdgeLen[1] + dualOppEdgeLen[0]*oppEdgeLen[0]);
    dualVertexVol[0] = 0.25 * (dualOppEdgeLen[1]*oppEdgeLen[1]*oppSin[1] + dualOppEdgeLen[2]*oppEdgeLen[2]*oppSin[2]);
    dualVertexVol[1] = 0.25 * (dualOppEdgeLen[0]*oppEdgeLen[0]*oppSin[0] + dualOppEdgeLen[2]*oppEdgeLen[2]*oppSin[2]);
    dualVertexVol[2] = 0.25 * (dualOppEdgeLen[1]*oppEdgeLen[1]*oppSin[1] + dualOppEdgeLen[0]*oppEdgeLen[0]*oppSin[0]);
  }


  //double ElVolumesInfo2d::getDiameter() {
  //  WorldVector<double> rVec = elInfo->getCoord(0) - cc;
  //  return 2.0 * sqrt(dot(rVec, rVec));
  //}
}
